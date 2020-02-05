import sys
import pickle
import h5py
import scipy as sp

if len(sys.argv) < 4:
    sys.stderr.write('Usage: %s <coordinate_db> <sparse_bam> <outfile> \n' % sys.argv[0])
    sys.exit(1)

junc_pickle = sys.argv[1]
sbam = sys.argv[2]
outfname = sys.argv[3]

STRANDS = ['none', 'p', 'm']

### load coordinate DB
junctions = pickle.load(open(junc_pickle, 'rb'))

### project sparse bam counts into space of  all junctions 
OUT = h5py.File(outfname, 'w')
junc_quant = dict()
with h5py.File(sbam, 'r') as IN:
    for key in IN:
        if key.endswith('_introns_m'):
            strand = 2
        elif key.endswith('_introns_p'):
            strand = 1
        else:
            continue
        chrm = key.split('_')[0]
        junc_quant = sp.zeros((len(junctions[chrm, strand]), ), dtype='int')

        #qidx = IN[key][:, 2] >= 3

        qidx = sp.where([(_[0], _[1]) in junctions[chrm, strand] for _ in IN[key][:]])[0]
        if sp.sum(qidx) == 0:
            outkey = 'junctions_%s_%s' % (chrm, STRANDS[strand])
            OUT.create_dataset(name=outkey, data=junc_quant[:, sp.newaxis], compression='gzip', maxshape=(junc_quant.shape[0], None)) 
            continue

        idxs = sp.array([junctions[chrm, strand][(_[0], _[1])] for _ in IN[key][qidx, :]])
        assert sp.all(idxs == sorted(idxs))
        junc_quant[idxs] =IN[key][qidx, 2] # curr_quant

        outkey = 'junctions_%s_%s' % (chrm, STRANDS[strand])
        OUT.create_dataset(name=outkey, data=junc_quant[:, sp.newaxis], compression='gzip', maxshape=(junc_quant.shape[0], None)) 
OUT.close()

