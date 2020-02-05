import sys
import os
import glob
import scipy as sp
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.markers as markers
from sklearn.decomposition import PCA

sys.path.append('../icgc_utils')
import names
# get lookup table of gene names
lookup = names.get_lookup_complete()

CONF = 2

sys.path.append('..')
from paths import BASEDIR,BASEDIR_AS

basedir = os.path.join(BASEDIR_AS, 'alternative_splicing')
metatable = os.path.join(BASEDIR, 'orig_data/metadata/per_aliquot_v2/rnaseq_metadata.histo.tsv')
plotdir = os.path.join(basedir, 'plots', 'psi_deviation_histo')
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

### prep metadata
metadata = []
for line in open(metatable, 'r'):
    sl = line.strip().split('\t')
    metadata.append(sl)
metadata = sp.array(metadata, dtype='str')
header = metadata[0, :]
metadata = metadata[1:, :]
pidx = sp.where(header == 'project_code')[0]
hidx = sp.where(header == 'histotype')[0]
alidx = sp.where(header == 'aliquot_id')[0]

infile_pat = ['merge_graphs_*_C%i.psi_dev_histo_noabs_per_gene_func.tsv.gz' % CONF, 'merge_graphs_*_C%i.psi_dev_histo_noabs_per_gene.tsv.gz' % CONF]
outfile_pat = ['merge_graphs_all_events_C%i.psi_dev_histo_noabs_per_gene_func.tsv.gz' % CONF, 'merge_graphs_all_events_C%i.psi_dev_histo_noabs_per_gene.tsv.gz' % CONF]
outfile_pat2 = ['merge_graphs_all_events_C%i.psi_dev_histo_noabs_per_gene_binary_func.tsv.gz' % CONF, 'merge_graphs_all_events_C%i.psi_dev_histo_noabs_per_gene_binary.tsv.gz' % CONF]
outfile_pat3 = ['merge_graphs_all_events_C%i.psi_dev_median_histo_noabs_per_gene_func.tsv.gz' % CONF, 'merge_graphs_all_events_C%i.psi_dev_median_histo_noabs_per_gene.tsv.gz' % CONF]
functional_tag = ['_func', '']

for p in [0]: #range(len(infile_pat)):

    filelist = glob.glob(os.path.join(basedir, infile_pat[p]))

    if 'plot_only' in sys.argv:
        print('Loading data from %s' % os.path.join(basedir, outfile_pat[p]))
        merged = sp.loadtxt(os.path.join(basedir, outfile_pat[p]), dtype='str', delimiter='\t')
        merged = merged[1:, :][:, 1:].astype('float')
        fname = filelist[0]
        print('processing %s for strains' % fname, file=sys.stderr)
        data = sp.loadtxt(fname, dtype='str', delimiter='\t')
        strains = data[0, 1:]
        del data
    else:
        if os.path.exists(os.path.join(basedir, outfile_pat[p])):
            os.remove(os.path.join(basedir, outfile_pat[p]))
        #filelist = glob.glob(os.path.join(basedir, infile_pat[p]))
        for f, fname in enumerate(filelist):
            print('processing %s' % fname, file=sys.stderr)
            data = sp.loadtxt(fname, dtype='str', delimiter='\t')
            data_ids = sp.loadtxt(re.sub(r'.tsv.gz$', '', fname) + '.event_ids.tsv.gz', dtype='str', delimiter='\t')
            if f == 0:
                genes = data[1:, 0]
                strains = data[0, 1:]
                merged = data[1:, :][:, 1:].astype('float')
                merged_ids = data_ids[1:, :][:, 1:]
                s_idx = sp.argsort(genes)
                genes = genes[s_idx]
                merged = merged[s_idx, :]
                merged_ids = merged_ids[s_idx, :]
                s_idx = sp.argsort(strains)
                strains = strains[s_idx]
                merged = merged[:, s_idx]
                merged_ids = merged_ids[:, s_idx]
            else:
                currgenes = data[1:, 0]
                currstrains = data[0, 1:]
                data = data[1:, :][:, 1:].astype('float')
                data_ids = data_ids[1:, :][:, 1:]
                s_idx = sp.argsort(currstrains)
                data = data[:, s_idx]
                data_ids = data_ids[:, s_idx]
                currstrains = currstrains[s_idx]
                assert sp.all(currstrains == strains)
                merged = sp.r_[merged, data]
                merged_ids = sp.r_[merged_ids, data_ids]
                genes = sp.r_[genes, currgenes]

                s_idx = sp.argsort(genes)
                genes = genes[s_idx]
                merged = merged[s_idx, :]
                merged_ids = merged_ids[s_idx, :]

                _, cnt = sp.unique(genes, return_counts=True)
                _, fidx = sp.unique(genes, return_index=True)
                newmerged = sp.empty((cnt.shape[0], merged.shape[1]), dtype='float')
                newmerged_ids = sp.empty((cnt.shape[0], merged_ids.shape[1]), dtype='|S40')
                cumcounts = 0
                for i,c in enumerate(cnt):
                    if c > 1:
                        tmp = sp.absolute(merged[cumcounts:cumcounts+c, :]).argmax(axis=0)
                        newmerged_ids[i, :] = sp.array([merged_ids[cumcounts:cumcounts+c, :][jj, ii] for ii, jj in enumerate(tmp)])
                        newmerged[i, :] = sp.array([merged[cumcounts:cumcounts+c, :][jj, ii] for ii, jj in enumerate(tmp)])
                    else:
                        newmerged[i, :] = merged[cumcounts, :]
                        newmerged_ids[i, :] = merged_ids[cumcounts, :]
                    cumcounts += c
                merged = newmerged
                merged_ids = newmerged_ids
                genes = genes[fidx]

        ### sort by mean deviation over strains
        s_idx = sp.argsort(sp.median(sp.absolute(merged), axis=1))[::-1]
        merged = merged[s_idx, :]
        merged_ids = merged_ids[s_idx, :]
        genes = genes[s_idx]

        print("Writing output files")
        merged_bin = (sp.absolute(merged) > 4.5).astype('int').astype('str')
        merged_bin = sp.r_[strains[sp.newaxis, :], merged_bin]
        merged_bin = sp.c_[sp.r_[['gene_id'], genes], merged_bin]

        merged_median = sp.median(sp.absolute(merged), axis=1).astype('str')
        merged_median = sp.append('median_dev', merged_median)
        merged_median = sp.c_[sp.append('gene_id', sp.array([x.split('.')[0] for x in genes])), merged_median]
        merged_median = sp.c_[sp.append('gene_name', sp.array([names.get_ID(x.split('.')[0], lookup) for x in genes])), merged_median]

        merged_str = merged.astype('str')
        merged_str = sp.r_[strains[sp.newaxis, :], merged_str]
        merged_str = sp.c_[sp.r_[['gene_id'], genes], merged_str]

        merged_ids_str = merged_ids.astype('str')
        merged_ids_str = sp.r_[strains[sp.newaxis, :], merged_ids_str]
        merged_ids_str = sp.c_[sp.r_[['gene_id'], genes], merged_ids_str]

        if not os.path.exists(os.path.join(basedir, outfile_pat[p])):
            sp.savetxt(os.path.join(basedir, outfile_pat[p]), merged_str, fmt='%s', delimiter='\t')
        if not os.path.exists(os.path.join(basedir, re.sub(r'.tsv.gz$', '', outfile_pat[p]) + '.event_ids.tsv.gz')):
            sp.savetxt(os.path.join(basedir, re.sub(r'.tsv.gz$', '', outfile_pat[p]) + '.event_ids.tsv.gz'), merged_ids_str, fmt='%s', delimiter='\t')
        #sp.savetxt(os.path.join(basedir, #outfile_pat2[p]), merged_bin, fmt='%s', delimiter='\t')
        if not os.path.exists(os.path.join(basedir, outfile_pat3[p])):
            sp.savetxt(os.path.join(basedir, outfile_pat3[p]), merged_median, fmt='%s', delimiter='\t')

    a,b = sp.where(strains[:, sp.newaxis] == metadata[:, alidx[0]])
    assert sp.all(strains == metadata[b, alidx[0]])
    projects = sp.array([':'.join([metadata[i, pidx[0]], metadata[i, hidx[0]]]) for i in b])
   # projects = metadata[b, pidx]
    projects_u = sp.unique(projects)

    print('Computing PCA')
    pca = PCA(n_components=2)
    trans_data = pca.fit_transform((merged.T > 4.5).astype('float'))
    cmap = plt.get_cmap('jet')
    norm = plt.Normalize(0, projects_u.shape[0])
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    for pp, project in enumerate(projects_u):
        curr_idx = sp.where(projects == project)[0]
        ax.plot(trans_data[curr_idx, 0], trans_data[curr_idx, 1], 'o', color=cmap(norm(pp)), alpha=0.75, label=project, marker=markers.MarkerStyle.filled_markers[pp % 13])
    #ax.plot(trans_data[:, 0], trans_data[:, 1], 'ob', alpha=0.75)
    ax.set_title('PCA - Integrated')
    ax.set_xlabel('PC 1')
    ax.set_ylabel('PC 2')
    ax.legend(numpoints=1, ncol=2, loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)

    plt.tight_layout()
    plt.savefig(os.path.join(plotdir, 'normalized_psi_pca_histo_noabs_integrated_C%i%s.png' % (CONF, functional_tag[p])), format='png', bbox_inches='tight')
    plt.close(fig)

