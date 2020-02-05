import sys
import os
import scipy as sp
import h5py
import cPickle

if len(sys.argv) < 2:
    print >> sys.stderr, 'Usage: %s <event_type>' % sys.argv[0]
    sys.exit(1)
event_type = sys.argv[1]

CONF = 2

sys.path.append('..')
from paths import BASEDIR_AS
basedir = os.path.join(BASEDIR_AS, 'alternative_splicing')

HGMD_DB = ### path to HGMD DB

if not os.path.exists('hgmd_pos.cpickle'):
    ### load database of HGMD mutations
    hgmd_pos = dict()
    for line in open(HGMD_DB, 'r'):
        if line[:2] == '##':
            continue
        sl = line.strip().split('\t')
        if line[0] == '#':
            ref_idx = sp.where(sp.array(sl, dtype='str') == 'REF')[0]
            alt_idx = sp.where(sp.array(sl, dtype='str') == 'ALT')[0]
            pos_idx = sp.where(sp.array(sl, dtype='str') == 'POS')[0]
            inf_idx = sp.where(sp.array(sl, dtype='str') == 'INFO')[0]
        else:
            ### only keep mutations that are classified as damaging
            tags = dict([x.split('=')[:2] for x in  sl[inf_idx].split(';')]) 
            if tags['CLASS'] != 'DM':
                continue

            varlen = max(len(sl[ref_idx]), len(sl[alt_idx]))
            try:
                hgmd_pos[sl[0]].extend(range(int(sl[pos_idx]), int(sl[pos_idx]) + varlen))
            except KeyError:
                hgmd_pos[sl[0]] = range(int(sl[pos_idx]), int(sl[pos_idx]) + varlen)
    #for i in range(hgmd.shape[0]):
    #    varlen = max(len(hgmd[i, ref_idx]), len(hgmd[i, alt_idx]))
    #    hgmd_pos.extend(range(int(hgmd[i, pos_idx]), int(hgmd[i, pos_idx]) + varlen))
    for k in hgmd_pos:
        hgmd_pos[k] = set(hgmd_pos[k])
    cPickle.dump(hgmd_pos, open('hgmd_pos.cpickle', 'w'), -1)
else:
    hgmd_pos = cPickle.load(open('hgmd_pos.cpickle', 'r'))


### load splicing events
IN = h5py.File(os.path.join(basedir, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, CONF)), 'r')
conf_idx = IN['conf_idx'][:].astype('int')

### walk through all events and flag the ones that overlap a HGMD variant
k1_idx = []
deltas = []
pos = IN['event_pos'][:]
for i, c in enumerate(conf_idx):
    if i > 0 and i % 1000 == 0:
        sys.stdout.write('.')
        if i % 10000 == 0:
            sys.stdout.write('%i/%i (%i)\n' % (i, conf_idx.shape[0], len(k1_idx)))
        sys.stdout.flush()
    chrm = IN['gene_chr'][IN['gene_idx'][c]]
    ### identify the delta positions of the event
    if event_type in ['alt_3prime', 'alt_5prime']:
        if pos[c, 4] == pos[c, 6]:
            delta = [pos[c, 1], pos[c, 3]]
        else:
            assert(pos[c, 1] == pos[c, 3])
            delta = [pos[c, 4], pos[c, 6]]
    elif event_type == 'exon_skip':
        delta = [pos[c, 2], pos[c, 3]]
    elif event_type == 'intron_retention':
        delta = [pos[c, 1], pos[c, 2]]
    elif event_type == 'mutex_exons':
        delta = [pos[c, 0, 1], pos[c, 1, 1], pos[c, 0, 2], pos[c, 1, 2]] 
    elif event_type == 'mult_exon_skip':
        delta = [pos[c, 2], pos[c, 3], pos[c, 4], pos[c, 5]] 
    deltas.append(delta)
    for d in range(len(delta) / 2):
        for j in range(delta[(2 * d) + 0], delta[(2 * d) + 1]):
            if j in hgmd_pos[chrm]:
                k1_idx.append(i)
                break
k1_idx = sp.unique(sp.array(k1_idx, dtype='int'))
print 'flagged %i of %i events that contain a variant annotated as disease causing by HGMD' % (k1_idx.shape[0], conf_idx.shape[0])

cPickle.dump(conf_idx[k1_idx], open(os.path.join(basedir, 'merge_graphs_%s_C%i.function_idx_hgmd_only.cpickle' % (event_type, CONF)), 'w'), -1)

### remove events that are divisible by three
k2_idx = []
assert len(deltas) == conf_idx.shape[0]
for i, c in enumerate(conf_idx):
    if i > 0 and i % 1000 == 0:
        sys.stdout.write('.')
        if i % 10000 == 0:
            sys.stdout.write('%i/%i\n' % (i, conf_idx.shape[0]))
        sys.stdout.flush()
    d = sp.array(deltas[i])
    #if (IN['event_pos'][c, 3] - IN['event_pos'][c, 2]) % 3 != 0:
    if sp.sum(d[1::2] - d[::2]) % 3 != 0:
        k2_idx.append(i)
k2_idx = sp.array(k2_idx, dtype='int')
print 'flagged %i of %i events that are out of frame' % (k2_idx.shape[0], conf_idx.shape[0])

### integrate k1 and k2 events
k_idx = sp.union1d(k1_idx, k2_idx)

conf_idx = conf_idx[k_idx]
#pos = pos[k_idx, :]
print 'retaining %i events' % (conf_idx.shape[0])

cPickle.dump(conf_idx, open(os.path.join(basedir, 'merge_graphs_%s_C%i.function_idx.cpickle' % (event_type, CONF)), 'w'), -1)

#print 'loading psi'
#psi = IN['psi'][:]
#print 'done'
#
#### remove all events that have nan as PSI in more than 10% of samples
#k_idx = []
#for i, c in enumerate(conf_idx):
#    if i > 0 and i % 1000 == 0:
#        sys.stdout.write('.')
#        if i % 10000 == 0:
#            sys.stdout.write('%i/%i\n' % (i, conf_idx.shape[0]))
#        sys.stdout.flush()
#    if sp.sum(sp.isnan(psi[:, c])).astype('float') / psi.shape[0] <= 0.1:
#        k_idx.append(i)
#k_idx = sp.array(k_idx, dtype='int')
#print 'remove %i of %i events due to high fraction of NaN' % (conf_idx.shape[0] - k_idx.shape[0], conf_idx.shape[0])
#conf_idx = conf_idx[k_idx]
#print 'retaining %i events' % (conf_idx.shape[0])

### remove events where no appreciable fraction of the samples is supported
#k_idx = []
#for i, c in enumerate(conf_idx):
#    if sp.sum(sp.absolute(psi[:, c] - 0.5) < 0.4) > (0.05 * psi.shape[0]):
#        k_idx.append(i)
#k_idx = sp.array(k_idx, dtype='int')
#print 'remove %i of %i events that have PSI between 0.05 and 0.95 in too few samples' % (conf_idx.shape[0] - k_idx.shape[0], conf_idx.shape[0])
#conf_idx = conf_idx[k_idx]
#pos = pos[k_idx, :]
#print 'retaining %i events' % (conf_idx.shape[0])
 
IN.close()
