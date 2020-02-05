import sys
import scipy as sp
import scipy.stats as spst
import h5py
import os
import re
import pickle

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.markers as markers
from sklearn.decomposition import PCA

CONF = 2

sys.path.append('..')
from paths import BASEDIR,BASEDIR_AS

basedir = os.path.join(BASEDIR_AS, 'alternative_splicing')
metatable = os.path.join(BASEDIR, 'ICGC/orig_data/metadata/per_aliquot_v2/rnaseq_metadata.histo.tsv')
plotdir = os.path.join(basedir, 'plots', 'psi_deviation_histo')
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

### settings
# min PSI deviation threshold
psi_thresh = 4
# do we compute the deviation from the PSI mean over all samples
# or relative to the mean of the individual subtypes?
per_subtype = True
# subset events to a set of functionally interesting ones
functional_only = False
functional_tag = ''
if functional_only:
    functional_tag = '_func'

event_short_dict = {'exon_skip':'es', 
                    'intron_retention':'ir',
                    'alt_3prime':'a3',
                    'alt_5prime':'a5',
                    'mult_exon_skip':'ms',
                    'mutex_exons':'mx'}

if len(sys.argv) < 2:
    print('Usage: %s <event_type>' % sys.argv[0], file=sys.stderr)
    sys.exit(1)
event_type = sys.argv[1]

### prep metadata
metadata = []
for line in open(metatable, 'r'):
    sl = line.strip().split('\t')
    metadata.append(sl)
metadata = sp.array(metadata, dtype='str')
header = metadata[0, :]
metadata = metadata[1:, :]
aidx = sp.where(header == 'analysis_id')[0]
alidx = sp.where(header == 'aliquot_id')[0]
tidx = sp.where(header == 'is_tumour')[0]
pidx = sp.where(header == 'project_code')[0]
hidx = sp.where(header == 'histotype')[0]

hdf5_countfile = os.path.join(basedir, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, CONF))
hdf5_psifile = os.path.join(basedir, 'merge_graphs_%s_C%i.counts.psi_conf.hdf5' % (event_type, CONF))
func_idx_file = os.path.join(basedir, 'merge_graphs_%s_C%i.function_idx.cpickle' % (event_type, CONF))

IN = h5py.File(hdf5_countfile, 'r')

### get relevant index vectors
cidx = IN['conf_idx'][:].astype('int')
gidx = IN['gene_idx'][:][cidx].astype('int')
strains = sp.array([x.split('.')[0] for x in IN['strains'][:]], dtype='str')
event_idx = cidx + 1

### load psi of confident events
if os.path.exists(hdf5_psifile):
    print('loading %i confident events from %s' % (cidx.shape[0], hdf5_psifile), file=sys.stderr)
    IN2 = h5py.File(hdf5_psifile, 'r')
    psi = IN2['psi'][:]
    IN2.close()
else:
    print('loading %i confident events from %s' % (cidx.shape[0], hdf5_countfile), file=sys.stderr)
    psi = sp.empty((IN['psi'].shape[0], cidx.shape[0]), dtype='float')
    for i in range(psi.shape[0]):
        if i > 0 and i % 10 == 0:
            sys.stderr.write('.')
            if i % 100 == 0:
                sys.stderr.write('%i/%i strains done\n' % (i, psi.shape[0]))
        tmp = IN['psi'][i, :]
        psi[i, :] = tmp[cidx]
    OUT = h5py.File(hdf5_psifile, 'w')
    OUT.create_dataset(name='psi', data=psi, compression='gzip')
    OUT.close()
print('\n...done', file=sys.stderr)

### sort metadata by strains and subset psi values to 
### donors present in metadata
a, b = sp.where(strains[:, sp.newaxis] == metadata[:, aidx[0]])
assert sp.all(strains[a] == metadata[:, aidx[0]][b]) 
strains = strains[a]
psi = psi[a, :]
metadata = metadata[b, :]

### remove normals from analysis
is_tumor = sp.where(metadata[:, tidx[0]] == 'yes')[0]
strains = strains[is_tumor]
psi = psi[is_tumor, :]
metadata = metadata[is_tumor, :]

### filter for functional events only
if functional_only:
    print('filtering for functional events', file=sys.stderr)
    cidx_filt = pickle.load(open(func_idx_file, 'r'))
    k_idx = sp.where(sp.in1d(cidx, cidx_filt))[0]
    psi = psi[:, k_idx]
    gidx = gidx[k_idx]
    cidx = cidx[k_idx]
    event_idx = event_idx[k_idx]
    print('retaining %s events' % k_idx.shape[0], file=sys.stderr)

### filter even more strictly
print('filtering for missing data', file=sys.stderr)
keep_mat = None
if per_subtype:
    #proj_idx = sp.where(header == 'histotype')[0][0]
    #projects = sp.unique(metadata[:, proj_idx])
    projects = sp.unique([':'.join([metadata[i, pidx[0]], metadata[i, hidx[0]]]) for i in range(metadata.shape[0])])
    keep_mat = sp.zeros((projects.shape[0], psi.shape[1]), dtype='bool')
    for p, project in enumerate(projects):
        print('processing subtype %s (%i/%i)' % (project, p + 1, projects.shape[0]))
        curr_idx = sp.where((metadata[:, pidx[0]] == project.split(':')[0]) & (metadata[:, hidx[0]] == project.split(':')[1]))[0]
        nan_frac = sp.mean(sp.isnan(psi[curr_idx, :]).astype('float'), axis=0)
        kidx = sp.where(nan_frac < 0.3)[0]
        keep_mat[p, kidx] = True
    kidx = keep_mat.max(axis=0)
else:
    #kidx = sp.where(sp.sum(sp.isnan(psi), axis=0) < (0.3 * psi.shape[0]))[0]
    kidx = sp.where(sp.mean(sp.isnan(psi).astype('float'), axis=0) < 0.3)[0]

psi = psi[:, kidx]
gidx = gidx[kidx]
cidx = cidx[kidx]
event_idx = event_idx[kidx]
if not keep_mat is None:    
    keep_mat = keep_mat[:, kidx]
print('%i events remain after filtering' % sp.sum(kidx), file=sys.stderr)
print('affecting %i genes' % sp.unique(gidx).shape[0], file=sys.stderr)
print('...done', file=sys.stderr)

### get dev from mean per subtype
if per_subtype:
    ### generate PSI distribution plots
    p_num = projects.shape[0]
    panels = sp.ceil(sp.sqrt(p_num)).astype('int')
    gs = gridspec.GridSpec(panels, panels)
    fig = plt.figure(figsize=(4*panels, 4*panels), dpi=200) 

    for p, project in enumerate(projects):
        print('computing mean for subtype %s (%i/%i)' % (project, p + 1, projects.shape[0]))
        ### operate on a part of the matrix
        #curr_idx = sp.where(metadata[:, proj_idx] == project)[0]
        curr_idx = sp.where((metadata[:, pidx[0]] == project.split(':')[0]) & (metadata[:, hidx[0]] == project.split(':')[1]))[0]
        kidx = sp.where(keep_mat[p, :])[0]
        kidxn = sp.where(~keep_mat[p, :])[0] 
        psi_mean = spst.nanmean(psi[curr_idx, :][:, kidx], axis=0)
        psi_stddev = spst.nanstd(psi[curr_idx, :][:, kidx], axis=0)
        for cc in curr_idx:
            psi[cc, kidx] = sp.absolute(psi[cc, kidx] - psi_mean)
            psi[cc, kidx] /= psi_stddev
            kkidx = sp.where(psi_stddev < 0.01)[0]
            psi[cc, kidx[kkidx]] = 0
            psi[cc, kidxn] = 0

        ax = fig.add_subplot(gs[p / panels, p % panels])
        print('plotting distribution for %s' % project)
        #ax.hexbin(psi[curr_idx, :][:, kidx].ravel(), sp.dot(sp.arange(min(curr_idx.shape[0], 100))[:, sp.newaxis], sp.ones((1, kidx.shape[0]))).ravel(), gridsize=(100, min(curr_idx.shape[0], 100)))
        #ax.hexbin(psi[curr_idx, :][:, kidx].ravel(), sp.dot(sp.arange(min(curr_idx.shape[0], 100))[:, sp.newaxis], sp.ones((1, kidx.shape[0]))).ravel(), gridsize=(10, min(curr_idx.shape[0], 100)), mincnt=1, bins='log')
        plt_mat = []
       # bounds = sp.array([ 0.01,  0.25,  0.5 ,  0.75,  1.  ,  1.25,  1.5 ,  1.75,  2.  ,
       #                     2.25,  2.5 ,  2.75,  3.  ,  3.25,  3.5 ,  3.75,  4.  ,  4.25,
       #                     4.5 ,  4.75,  5.  ], dtype='float')
        bounds = sp.array([ 0.01,  0.1 ,  0.25,  0.5 ,  1.,  2.,  4.,  8.,  16 ], dtype='float')

        for cc in curr_idx:
            nnidx = sp.where(~sp.isnan(psi[cc, kidx]))[0]
            bins, _ = sp.histogram(psi[cc, kidx][nnidx], bins=bounds)
            ax.plot(sp.arange(bins.shape[0]) + 0.5, bins, 'b-')
            #plt_mat.append(bins)
        #plt_mat = sp.array(plt_mat)
        #ax.imshow(sp.sqrt(plt_mat), aspect='auto', interpolation='none')
        #prct95 = spst.scoreatpercentile(psi[curr_idx, :][:, kidx].ravel(), 95)
        #ylim = ax.get_ylim()
        #xlim = ax.get_xlim()
        #ax.plot([prct95, prct95], ylim, '--r', linewidth=2)
        #ax.set_ylim(ylim)
        #ax.set_xlim(xlim)
        ax.set_title(project + '(N=%i)' % curr_idx.shape)
        ax.set_xticks(sp.arange(bounds.shape[0]))
        ax.set_xticklabels(bounds, rotation=90)
        #ax.set_xticks(sp.arange(20)[::2] + 0.5)
        #ax.set_xticklabels(bounds[::2][:-1] + 0.25, rotation=90)
        ax.set_xlabel('PSI deviation')
        ax.set_ylabel('Frequency')

    plt.tight_layout()
    plt.savefig(os.path.join(plotdir, 'normalized_psi_distribution_histo_%s%s.png' % (event_type, functional_tag)), format='png', bbox_inches='tight')
    plt.close(fig)
        

### get dev from mean over all samples
else:
    psi_mean = spst.nanmean(psi, axis=0)
    psi_stddev = spst.nanstd(psi, axis=0)
    psi = sp.absolute(psi - sp.tile(psi_mean, (psi.shape[0], 1)))
    psi /= sp.tile(psi_stddev, (psi.shape[0], 1))

### substitute NaN values
psi[sp.isnan(psi)] = 0

#import pdb
#pdb.set_trace()

### compute PCA on psi table
print('Computing Event Level PCA')
pca = PCA(n_components=2)
trans_data = pca.fit_transform(psi)
cmap = plt.get_cmap('jet')
norm = plt.Normalize(0, projects.shape[0])
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
for p, project in enumerate(projects):
    #curr_idx = sp.where(metadata[:, proj_idx] == project)[0]
    curr_idx = sp.where((metadata[:, pidx[0]] == project.split(':')[0]) & (metadata[:, hidx[0]] == project.split(':')[1]))[0]
    ax.plot(trans_data[curr_idx, 0], trans_data[curr_idx, 1], 'o', color=cmap(norm(p)), alpha=0.75, label=project, marker=markers.MarkerStyle.filled_markers[p % 13])
    ax.set_title('PCA - Event Level (%s)' % event_type)
    ax.set_xlabel('PC 1')
    ax.set_ylabel('PC 2')
ax.legend(numpoints=1, ncol=2, loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)

plt.tight_layout()
plt.savefig(os.path.join(plotdir, 'normalized_psi_pca_histo_%s%s.png' % (event_type, functional_tag)), format='png', bbox_inches='tight')
plt.close(fig)

### collapse per gene
### - if there are multiple events per gene, take event with max
###   PSI to represent the full gene
gsidx = sp.argsort(gidx)
gidx = gidx[gsidx]
gidx,fidx = sp.unique(gidx, return_index=True)
lidx = sp.r_[fidx[1:], gsidx.shape[0]]
gene_names = IN['gene_names'][:][gidx]

psi_collapsed = sp.empty((psi.shape[0], fidx.shape[0]), dtype='float')
event_idx_collapsed = sp.empty((psi.shape[0], fidx.shape[0]), dtype='int')
for i in range(lidx.shape[0]):
    psi_collapsed[:, i] = psi[:, gsidx[fidx[i]:lidx[i]]].max(axis=1)
    tmp = psi[:, gsidx[fidx[i]:lidx[i]]].argmax(axis=1)
    event_idx_collapsed[:, i] = event_idx[gsidx[fidx[i]:lidx[i]]][tmp]

### gene centric PCA
print('Computing Gene Level PCA')
pca = PCA(n_components=2)
trans_data = pca.fit_transform(psi_collapsed)
cmap = plt.get_cmap('jet')
norm = plt.Normalize(0, projects.shape[0])
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
for p, project in enumerate(projects):
    #curr_idx = sp.where(metadata[:, proj_idx] == project)[0]
    curr_idx = sp.where((metadata[:, pidx[0]] == project.split(':')[0]) & (metadata[:, hidx[0]] == project.split(':')[1]))[0]
    ax.plot(trans_data[curr_idx, 0], trans_data[curr_idx, 1], 'o', color=cmap(norm(p)), alpha=0.75, label=project, marker=markers.MarkerStyle.filled_markers[p % 13])
    ax.set_title('PCA - Gene Level (%s)' % event_type)
    ax.set_xlabel('PC 1')
    ax.set_ylabel('PC 2')
ax.legend(numpoints=1, ncol=2, loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)

plt.tight_layout()
plt.savefig(os.path.join(plotdir, 'normalized_psi_per_gene_pca_histo_%s%s.png' % (event_type, functional_tag)), format='png', bbox_inches='tight')
plt.close(fig)

psi_collapsed = psi_collapsed.T
psi_collapsed_bin = (psi_collapsed >= psi_thresh)
event_idx_collapsed = event_idx_collapsed.T.astype('str')
for i in range(event_idx_collapsed.shape[0]):
    for j in range(event_idx_collapsed.shape[1]):
        event_idx_collapsed[i, j] = event_short_dict[event_type] + ':' + event_idx_collapsed[i, j]

psi_collapsed = psi_collapsed.astype('str')
#psi_collapsed = sp.r_[strains[sp.newaxis, :], psi_collapsed]
psi_collapsed = sp.r_[metadata[:, alidx[0]][sp.newaxis, :], psi_collapsed]
psi_collapsed = sp.c_[sp.hstack(['gene_id', gene_names]), psi_collapsed]
outfname = os.path.join(basedir, 'merge_graphs_%s_C%i.psi_dev_histo_per_gene%s.tsv.gz' % (event_type, CONF, functional_tag))
if not os.path.exists(outfname):
    sp.savetxt(outfname, psi_collapsed, fmt='%s', delimiter='\t')

event_idx_collapsed = sp.r_[metadata[:, alidx[0]][sp.newaxis, :], event_idx_collapsed]
event_idx_collapsed = sp.c_[sp.hstack(['gene_id', gene_names]), event_idx_collapsed]
outfname = os.path.join(basedir, 'merge_graphs_%s_C%i.psi_dev_histo_per_gene%s.event_ids.tsv.gz' % (event_type, CONF, functional_tag))
if not os.path.exists(outfname):
    print('writing ' + outfname)
    sp.savetxt(outfname, event_idx_collapsed, fmt='%s', delimiter='\t')

kidx = sp.where(sp.sum(psi_collapsed_bin, axis=1) > 0)[0] 
psi_collapsed_bin = psi_collapsed_bin[kidx, :]
gene_names = gene_names[kidx]

### how sparse is the binary matrix for the given cutoff?
sparsity = sp.sum(psi_collapsed_bin.astype('int') != 0).astype('float') / (psi_collapsed_bin.shape[0] * psi_collapsed_bin.shape[1]) * 100
print('given the current threshold of %i, %.2f percent of the matrix are non-zero' % (psi_thresh, sparsity))

psi_collapsed_bin = psi_collapsed_bin.astype('int').astype('str')
#psi_collapsed_bin = sp.r_[strains[sp.newaxis, :], psi_collapsed_bin]
psi_collapsed_bin = sp.r_[metadata[:, alidx[0]][sp.newaxis, :], psi_collapsed_bin]
psi_collapsed_bin = sp.c_[sp.hstack(['gene_id', gene_names]), psi_collapsed_bin]
outfname = os.path.join(basedir, 'merge_graphs_%s_C%i.psi_dev_histo_per_gene_binary%s.min%i.tsv.gz' % (event_type, CONF, functional_tag, psi_thresh))
sp.savetxt(outfname, psi_collapsed_bin, fmt='%s', delimiter='\t')

IN.close()
