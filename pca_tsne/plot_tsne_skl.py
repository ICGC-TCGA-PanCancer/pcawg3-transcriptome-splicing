import scipy as sp
import scipy.stats as spst
from scipy.linalg import eigh
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
import os
import pdb
import utils
import matplotlib.markers as markers
import pickle
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

sys.path.append('/cluster/home/akahles/git/tools/python')
from viz.axes import *

sys.path.append('/cluster/home/akahles/git/projects/2014/icgc/colors')
import ICGC_colors as ic

from paths import BASEDIR,BASEDIR_AS

event_label_dict = {'exon_skip' : 'Exon Skip',
                    'intron_retention' : 'Intron Retention',
                    'alt_3prime' : 'Alternative 3\'',
                    'alt_5prime' : 'Alternative 5\'',
                    'mutex_exons' : 'Mutually Exclusive Exons',
                    'mult_exon_skip' : 'Coordinated Exon Skip'}

event_types = ['exon_skip', 'alt_5prime', 'alt_3prime', 'intron_retention', 'mutex_exons', 'mult_exon_skip']
basedir_icgc = os.path.join(BASEDIR_AS, 'alternative_splicing')
basedir_gtex = os.path.join(BASEDIR_AS, 'alternative_splicing_gtex_count_only')
metadata_icgc  = os.path.join(BASEDIR, 'orig_data/metadata/per_aliquot_v2/rnaseq_metadata.tsv')
histodata_icgc = os.path.join(BASEDIR, 'qtl_analysis/annotation/pcawg_specimen_histology_August2016_v6_edited.tsv')
coding_genes = sp.loadtxt(os.path.join(BASEDIR, 'annotation/gencode.v19.annotation.hs37d5_chr.gtf.coding_genes.txt'), delimiter='\t', dtype='str')
expression_tsv = os.path.join(BASEDIR, 'librarySizeNorm/gene_counts_UCSC_freeze_7_aliquotlevel/graylist/joint_fpkm_uq.tsv.gz')
expression_filter_pickle = 'gene_exp_no_outliers.pickle'

datadir = os.path.join(basedir_icgc, 'pca')
if not os.path.exists(datadir):
    os.makedirs(datadir)
plotdir = os.path.join(basedir_icgc, 'plots', 'tsne')
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

if '-h' in sys.argv or '--help' in sys.argv:
    print('Usage: %s [<conf> (0.01)]' % sys.argv[0], file=sys.stderr)
    sys.exit(1)
k = 2
perplexity = 100

if len(sys.argv) < 2:
    conf = 0.01
elif sys.argv[1][0] != '-':
    conf = float(sys.argv[1])
else:
    conf = 0.01

filter_tag = ''
do_filter = ('--filter' in sys.argv)
if do_filter:
    filter_tag = '.filtered'

# create ICGC cancer type dictionary
(ct_dict, is_tumor_dict) = utils.get_ct_dict_metatable(metadata_icgc)

# create ICGC histotype dictionary
(ht_dict) = utils.get_ht_dict_metatable(metadata_icgc, histodata_icgc)

# create GTex type dictionary
sample_dict = {'Bladder-GTEx' :os.path.join(BASEDIR, 'gtex_tables/GTex_Bladder_samples.txt'),
               'Breast-GTEx'  :os.path.join(BASEDIR, 'gtex_tables/GTex_Breast_samples.txt'),
               'Kidney-GTEx'  :os.path.join(BASEDIR, 'gtex_tables/GTex_Kidney_samples.txt'),
               'Liver-GTEx'   :os.path.join(BASEDIR, 'gtex_tables/GTex_Liver_samples.txt'),
               'Lung-GTEx'    :os.path.join(BASEDIR, 'gtex_tables/GTex_Lung_samples.txt'),
               'Prostate-GTEx':os.path.join(BASEDIR, 'gtex_tables/GTex_Prostate_samples.txt'),
               'Stomache-GTEx':os.path.join(BASEDIR, 'gtex_tables/GTex_Stomach_samples.txt'),
               'Thyroid-GTEx' :os.path.join(BASEDIR, 'gtex_tables/GTex_Thyroid_samples.txt'),
               'Uterus-GTEx'  :os.path.join(BASEDIR, 'gtex_tables/GTex_Uterus_samples.txt'),
               'REM.'         :os.path.join(BASEDIR, 'gtex_tables/GTex_rest_wo_cells_samples.txt')}
gt_dict = utils.get_gt_dict(sample_dict)

if not os.path.exists(expression_filter_pickle):
    print('loading expression data from %s' % expression_tsv, file=sys.stderr)
    exp_data = sp.loadtxt(expression_tsv, delimiter='\t', dtype='str')
    exp_genes = exp_data[1:, :][:, 0]
    exp_data = exp_data[1:, :][:, 1:].astype('float') 
    exp_med = sp.median(exp_data, axis=1)
    k_idx = sp.where(exp_med > spst.scoreatpercentile(exp_med, 50))[0]
    exp_genes = exp_genes[k_idx]
    pickle.dump(exp_genes, open(expression_filter_pickle, 'w'), -1)
else:
    exp_genes = pickle.load(open(expression_filter_pickle, 'r'))

for event_type in event_types[:1]:
    print('processing %s' % event_type)

    picklefile = '%s/tsne_skl_%s.TN.conf_%.2f.perp_%i%s.pickle' % (datadir, event_type, 1.0 - conf, perplexity, filter_tag)

    if not os.path.exists(picklefile):

        ### get indices of confident events 
        IN = h5py.File('%s/merge_graphs_%s_C3.counts.hdf5' % (basedir_icgc, event_type), 'r')
        c_idx = IN['conf_idx'][:].astype('int')
        IN.close()
        IN = h5py.File('%s/merge_graphs_%s_C3.counts.hdf5' % (basedir_gtex, event_type), 'r')
        c_idx_gt = IN['conf_idx'][:].astype('int') 
        c_idx_gt = sp.intersect1d(c_idx_gt, c_idx)
        c_idx = sp.intersect1d(c_idx, c_idx_gt)
        IN.close()
        assert sp.all(c_idx == c_idx_gt)

        ### load ICGC data from hdf5
        print('Loading data from ICGC hdf5')
        IN = h5py.File('%s/merge_graphs_%s_C3.counts.hdf5' % (basedir_icgc, event_type), 'r')
        c_idx = IN['conf_idx'][:].astype('int')
        strains = sp.array([x.split('.')[0] for x in IN['strains'][:]], dtype='str')

        ### get psi values
        psi = sp.empty((IN['psi'].shape[0], c_idx.shape[0]), dtype='float')
        chunksize = IN['psi'].chunks[1] * 30
        cum = 0
        for c, chunk in enumerate(range(0, IN['psi'].shape[1], chunksize)):
            sys.stdout.write('.')
            if c > 0 and c % 10 == 0:
                sys.stdout.write('%i/%i chunks done\n' % (c, IN['psi'].shape[1] / chunksize + 1))
            sys.stdout.flush()
            curr_range = list(range(chunk, min(chunk + chunksize, IN['psi'].shape[1])))
            k_idx = sp.where(sp.in1d(curr_range, c_idx))[0]
            if k_idx.shape[0] > 0:
                tmp = IN['psi'][:, chunk:min(chunk + chunksize, IN['psi'].shape[1])]
                psi[:, cum:(cum + k_idx.shape[0])] = tmp[:, k_idx]
                cum += k_idx.shape[0]

        ### only keep events in coding genes
        gene_names = IN['gene_names'][:]
        gene_idx = IN['gene_idx'][:].astype('int')
        genes = sp.array([gene_names[gene_idx[_]] for _ in IN['conf_idx'][:]]) 
        coding_idx = sp.where(sp.in1d(genes, coding_genes))[0]
        IN.close()

        ### filter to coding genes
        psi = psi[:, coding_idx]
        genes = genes[coding_idx]

        if do_filter:
            ### only keep events in highly expressed genes
            exp_idx = sp.where(sp.in1d(genes, exp_genes))[0]
            psi = psi[:, exp_idx]
            genes = genes[exp_idx]

            ### only keep events in genes that are not in the top 5% of
            ### most complicated splice graph genes (measured on number of edges)
            genes_u, gene_counts = sp.unique(genes, return_counts=True)
            comp_idx = sp.where(gene_counts < spst.scoreatpercentile(gene_counts, 95))[0]
            comp_idx = sp.where(sp.in1d(genes, genes_u[comp_idx]))[0]
            psi = psi[:, comp_idx]
            genes = genes[comp_idx]

        ### only keep strains that are present in the metadata
        k_idx = sp.where(sp.in1d(strains, list(ct_dict.keys())))[0]
        strains = strains[k_idx]
        psi = psi[k_idx, :]

        ### sort by strains
        s_idx = sp.argsort(strains)
        strains = strains[s_idx]
        psi = psi[s_idx, :]

        ### subset to samples where we have normals
        #ct_idx = sp.where(sp.in1d(ctypes, ctypes[sp.where(tn_labels == 'Normal')[0]]))[0]
        #strains = strains[ct_idx]
        #psi = psi[ct_idx, :]


        ### load GTEx data from hdf5
        print('Loading data from GTEx hdf5')
        IN = h5py.File('%s/merge_graphs_%s_C3.counts.hdf5' % (basedir_gtex, event_type), 'r')
        c_idx_gt = IN['conf_idx'][:].astype('int') 
        c_idx_gt = sp.intersect1d(c_idx_gt, c_idx)

        strains_gt = sp.array([x.split('.')[0] for x in IN['strains'][:]])
        ### get psi values
        psi_gt = sp.empty((IN['psi'].shape[0], c_idx.shape[0]), dtype='float')
        chunksize = IN['psi'].chunks[1] * 30
        cum = 0
        for c, chunk in enumerate(range(0, IN['psi'].shape[1], chunksize)):
            sys.stdout.write('.')
            if c > 0 and c % 10 == 0:
                sys.stdout.write('%i/%i chunks done\n' % (c, IN['psi'].shape[1] / chunksize + 1))
            sys.stdout.flush()
            curr_range = list(range(chunk, min(chunk + chunksize, IN['psi'].shape[1])))
            k_idx = sp.where(sp.in1d(curr_range, c_idx))[0]
            if k_idx.shape[0] > 0:
                tmp = IN['psi'][:, chunk:min(chunk + chunksize, IN['psi'].shape[1])]
                psi_gt[:, cum:(cum + k_idx.shape[0])] = tmp[:, k_idx]
                cum += k_idx.shape[0]
        IN.close()

        ### remove everything we don't want
        k_idx = sp.where(sp.in1d(strains_gt, list(gt_dict.keys())))[0]
        strains_gt = strains_gt[k_idx]
        psi_gt = psi_gt[k_idx]

        ### filter to coding genes
        psi_gt = psi_gt[:, coding_idx]

        if do_filter:
            ### filter to highly expressed genes
            psi_gt = psi_gt[:, exp_idx]

            ### filter to highly complex genes
            psi_gt = psi_gt[:, comp_idx]

        ### remove everything from category REM
        ctypes_gt = sp.array([gt_dict[x] for x in strains_gt])
        k_idx = sp.where(ctypes_gt != 'REM.')[0]
        strains_gt = strains_gt[k_idx]
        psi_gt = psi_gt[k_idx, :]

        psi = sp.r_[psi, psi_gt].T

        ### filter for NaN
        k_idx = sp.where(sp.sum(sp.isnan(psi), axis=1) < conf * psi.shape[1])[0]

        psi = psi[k_idx, :]
        print('keeping %i events' % k_idx.shape[0])

        ### mean imputation for remaining nans
        print('Imputing mean')
        for i in range(psi.shape[0]):
            n_idx = sp.where(sp.isnan(psi[i, :]))[0]
            if n_idx.shape[0] == 0:
                continue
            psi[i, n_idx] = spst.nanmean(psi[i, :])

        ### center the data - I might not need to do this for the covariance
        psi -= sp.mean(psi, axis=1)[:, sp.newaxis]

        ### compute kernel
        #print 'Computing covariances'
        #K = sp.cov([psi[i, :] for i in range(psi.shape[0])])

        ### PCA 
        print('Compute t-SNE ...')
        #w_g, Vt_g = eigh(K)
        #V_g = Vt_g.T
        #w_g = w_g[::-1]
        #V_g = V_g[::-1, :]
        tsne = TSNE(n_components=2, init='pca', random_state=0, perplexity=perplexity)
        trans_data = tsne.fit_transform(psi.T)
        print('... done')

        #cPickle.dump((w_g, V_g, ctypes, tn_labels, psi), open(picklefile, 'w'), -1)
        pickle.dump((trans_data, strains, strains_gt), open(picklefile, 'w'), -1)
    else:
        print('Loading data from pickle: %s' % picklefile)
        #(w_g, V_g, ctypes, tn_labels, psi) = cPickle.load(open(picklefile, 'r'))
        (trans_data, strains, strains_gt) = pickle.load(open(picklefile, 'r'))


    ### get sample labels
    ctypes = sp.array([ct_dict[x] for x in strains])
    tn_labels = sp.array(['Tumor' if is_tumor_dict[x] else 'Normal' for x in strains])
    htypes = sp.array([ht_dict[x] for x in strains])
    nn_idx = (htypes!='NA')
    htype_colors = dict(list(zip(sp.unique(htypes[nn_idx]), ic.get_color_scheme('tumor_subtype', labels=sp.unique(htypes[nn_idx])))))
    ctypes_gt = sp.array([gt_dict[x] for x in strains_gt])
    tn_labels = sp.r_[tn_labels, sp.array(['GTEx' for x in strains_gt])]
    strains = sp.r_[strains, strains_gt]
    ctypes = sp.r_[ctypes, ctypes_gt]
    htypes_u = sp.unique(htypes)
    ctypes_gt_u = sp.unique(ctypes_gt)
    htypes = sp.r_[htypes, ctypes_gt]

    trans_data = trans_data.T

    ### choose colormap and adapt normalization 
    cmap = plt.get_cmap('jet')
    norm = plt.Normalize(0, sp.unique(ctypes).shape[0])
    ### plot first k main axes of variation
    print('Plotting by ctype ... ')
    fig = plt.figure(figsize = (k*4, k*4))
    for k1 in range(0, k):
        cnt = 1
        for k2 in range(k1 + 1, k):
            ax = fig.add_subplot(k-1, k-1, (k1 * (k-1)) + cnt)
            cnt += 1
            #trans_data = sp.vstack([V_g[k1, :], V_g[k2, :]]).dot(psi)
            for idx, ct in enumerate(sp.unique(ctypes)):
                c_idx = sp.where(ctypes == ct)[0]
                if c_idx.shape[0] > 0:
                    ax.plot(trans_data[k1, c_idx], trans_data[k2, c_idx], markers.MarkerStyle.filled_markers[idx % 13], color = cmap(norm(idx)), label = ct, ms = 8, alpha=0.75)
            ax.set_title('t-SNE Visualization for %s Events' % event_label_dict[event_type])
            ax.set_xticks(ax.get_xticks()[::2])
            ax.set_yticks(ax.get_yticks()[::2])
            ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
            ax.tick_params(axis = 'both', which = 'minor', labelsize = 10)

            clean_axis(ax)
            set_ticks_outer(ax)
        #if k1 == (k - 1):
        #ax.legend(numpoints = 1, ncol = 2, loc = 'center left', bbox_to_anchor = (1.2, 0.5))
        legend = ax.legend(numpoints = 1, ncol = 2, loc = 'center left', bbox_to_anchor = (1.1, 0.6), frameon=False)
        for l in legend.legendHandles:
            l._markersize = 10


    plt.tight_layout()
    plt.savefig(os.path.join(plotdir, 'tsne_skl_%s.CT.%i.conf_%.2f.perp_%i%s.png' % (event_type, k, 1.0 - conf, perplexity, filter_tag)), dpi = 300, format = 'png', bbox_inches='tight')
    plt.savefig(os.path.join(plotdir, 'tsne_skl_%s.CT.%i.conf_%.2f.perp_%i%s.pdf' % (event_type, k, 1.0 - conf, perplexity, filter_tag)), dpi = 300, format = 'pdf', bbox_inches='tight')
    plt.close(fig)
    print('... done.')

    ### choose colormap and adapt normalization 
    cmap = plt.get_cmap('jet')
    #norm = plt.Normalize(0, sp.unique(htypes).shape[0])
    norm = plt.Normalize(0, ctypes_gt_u.shape[0])
    ### plot first k main axes of variation
    print('Plotting by htype ... ')
    fig = plt.figure(figsize = (k*4, k*4))
    for k1 in range(0, k):
        cnt = 1
        for k2 in range(k1 + 1, k):
            ax = fig.add_subplot(k-1, k-1, (k1 * (k-1)) + cnt)
            cnt += 1
            #trans_data = sp.vstack([V_g[k1, :], V_g[k2, :]]).dot(psi)
            for idx, ct in enumerate(htypes_u):
                if ct == 'NA':
                    continue
                c_idx = sp.where(htypes == ct)[0]
                try:
                    color = htype_colors[ct]
                except KeyError:
                    color = cmap(norm(idx))
                if c_idx.shape[0] > 0:
                    ax.plot(trans_data[k1, c_idx], trans_data[k2, c_idx], markers.MarkerStyle.filled_markers[idx % 13], color=color, label=ct, ms=8, alpha=0.75)
            for idx, ct in enumerate(ctypes_gt_u):
                if ct == 'NA':
                    continue
                c_idx = sp.where(htypes == ct)[0]
                try:
                    color = htype_colors[ct]
                except KeyError:
                    color = cmap(norm(idx))
                if c_idx.shape[0] > 0:
                    ax.plot(trans_data[k1, c_idx], trans_data[k2, c_idx], markers.MarkerStyle.filled_markers[0], color=color, label=ct, ms=8, alpha=0.75)
            c_idx = sp.where(htypes == 'NA')[0]
            ax.plot(trans_data[k1, c_idx], trans_data[k2, c_idx], markers.MarkerStyle.filled_markers[sp.where(htypes_u == 'NA')[0][0] % 13], color = '0.6', label = 'Normal tissue', ms = 8, alpha=0.5)
            ax.set_title('t-SNE Visualization for %s Events' % event_label_dict[event_type])
            ax.set_xticks(ax.get_xticks()[::2])
            ax.set_yticks(ax.get_yticks()[::2])
            ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
            ax.tick_params(axis = 'both', which = 'minor', labelsize = 10)

            clean_axis(ax)
            set_ticks_outer(ax)
        #if k1 == (k - 1):
        #ax.legend(numpoints = 1, ncol = 2, loc = 'center left', bbox_to_anchor = (1.2, 0.5))
        legend = ax.legend(numpoints = 1, ncol = 2, loc = 'center left', bbox_to_anchor = (1.1, 0.6), frameon=False)
        for l in legend.legendHandles:
            l._markersize = 10

    plt.tight_layout()
    plt.savefig(os.path.join(plotdir, 'tsne_skl_%s.HT.%i.conf_%.2f.perp_%i%s.png' % (event_type, k, 1.0 - conf, perplexity, filter_tag)), dpi = 300, format = 'png', bbox_inches='tight')
    plt.savefig(os.path.join(plotdir, 'tsne_skl_%s.HT.%i.conf_%.2f.perp_%i%s.pdf' % (event_type, k, 1.0 - conf, perplexity, filter_tag)), dpi = 300, format = 'pdf', bbox_inches='tight')
    plt.close(fig)
    print('... done.')


    ### choose colormap and adapt normalization 
    cmap = plt.get_cmap('jet')
    norm = plt.Normalize(0, sp.unique(tn_labels).shape[0])
    ### plot first k main axes of variation
    print('Plotting by tumor / normal ... ')
    fig = plt.figure(figsize = (k*4, k*4))
    for k1 in range(0, k):
        cnt = 1
        for k2 in range(k1 + 1, k):
            ax = fig.add_subplot(k-1, k-1, (k1 * (k-1)) + cnt)
            cnt += 1
            #trans_data = sp.dot(sp.vstack([V_g[k1, :], V_g[k2, :]]), psi)
            for idx, ct in enumerate(sp.unique(tn_labels)[::-1]):
                c_idx = sp.where(tn_labels == ct)[0]
                if c_idx.shape[0] > 0:
                    #ax.plot(V_g[k1, c_idx], V_g[k2, c_idx], 'o', color = cmap(norm(idx)), label = ct, ms = 4, alpha=0.75)
                    ax.plot(trans_data[k1, c_idx], trans_data[k2, c_idx], 'o', color = cmap(norm(idx)), label = ct, ms = 10, alpha=0.75)
            ax.set_title('t-SNE Visualization for %s Events' % event_label_dict[event_type])
            ax.set_xticks(ax.get_xticks()[::2])
            ax.set_yticks(ax.get_yticks()[::2])
            ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
            ax.tick_params(axis = 'both', which = 'minor', labelsize = 10)

            clean_axis(ax)
            set_ticks_outer(ax)
        #if k1 == (k - 1):
        legend = ax.legend(numpoints = 1, ncol = 2, loc = 'center left', bbox_to_anchor = (1.1, 0.6), frameon=False)

    plt.tight_layout()

    plt.savefig(os.path.join(plotdir, 'tsne_skl_%s.TN.%i.conf_%.2f.perp_%i%s.png' % (event_type, k, 1.0 - conf, perplexity, filter_tag)), dpi = 300, format = 'png', bbox_inches='tight')
    plt.savefig(os.path.join(plotdir, 'tsne_skl_%s.TN.%i.conf_%.2f.perp_%i%s.pdf' % (event_type, k, 1.0 - conf, perplexity, filter_tag)), dpi = 300, format = 'pdf', bbox_inches='tight')
    print('... done.')
    plt.close(fig)

    print('Plotting individual ctype plots ...')
    sets = ['bladder', 'breast', 'lung', 'kidney', 'prostate', 'thyroid', 'stomache', 'uterus', 'liver']
    for s in sets:
        print('... %s' % s)
        fig = plt.figure(figsize = (k*4, k*4))
        for k1 in range(0, k):
            cnt = 1
            for k2 in range(k1 + 1, k):
                ax = fig.add_subplot(k-1, k-1, (k1 * (k-1)) + cnt)
                cnt += 1
                #trans_data = sp.vstack([V_g[k1, :], V_g[k2, :]]).dot(psi) 
                ### plot active elements
                if s == 'bladder':
                    norm = plt.Normalize(0, 6)
                    ct1 = ['BLCA-US']
                    ct2 = ['Bladder-GTEx']
                elif s == 'breast':
                    norm = plt.Normalize(0, 6)
                    ct1 = ['BRCA-US']
                    ct2 = ['Breast-GTEx']
                elif s == 'kidney':
                    norm = plt.Normalize(0, 6)
                    ct1 = ['KICH-US', 'KIRC-US', 'KIRP-US']
                    ct2 = ['Kidney-GTEx']
                elif s == 'lung':
                    norm = plt.Normalize(0, 6)
                    ct1 = ['LUSC-US', 'LUAD-US']
                    ct2 = ['Lung-GTEx']
                elif s == 'thyroid':
                    norm = plt.Normalize(0, 6)
                    ct1 = ['THCA-US']
                    ct2 = ['Thyroid-GTEx']
                elif s == 'uterus':
                    norm = plt.Normalize(0, 6)
                    ct1 = ['UCEC-US']
                    ct2 = ['Uterus-GTEx']
                elif s == 'stomache':
                    norm = plt.Normalize(0, 6)
                    ct1 = ['STAD-US']
                    ct2 = ['Stomache-GTEx']
                elif s == 'prostate':
                    norm = plt.Normalize(0, 6)
                    ct1 = ['PRAD-US']
                    ct2 = ['Prostate-GTEx']
                elif s == 'liver':
                    norm = plt.Normalize(0, 6)
                    ct1 = ['LIHC-US', 'LIRI-JP']
                    ct2 = ['Liver-GTEx']
                c1_idx = sp.where(sp.array([_ in ct1 for _ in ctypes]) & (tn_labels == 'Tumor'))[0]
                c2_idx = sp.where(sp.array([_ in ct1 for _ in ctypes]) & (tn_labels == 'Normal'))[0]
                c3_idx = sp.where(sp.array([_ in ct2 for _ in ctypes]))[0]
                cc_idx = sp.hstack([c1_idx, c2_idx, c3_idx])
                c6_idx = sp.where(~sp.in1d(sp.arange(ctypes.shape[0]), cc_idx))[0]
                ### plot inactive elements
                if c6_idx.shape[0] > 0:
                    ax.plot(trans_data[k1, c6_idx], trans_data[k2, c6_idx], markers.MarkerStyle.filled_markers[5], color = '0.6', label='other', ms = 10, alpha=0.5)
                ### plot active elements
                if c1_idx.shape[0] > 0:
                    ax.plot(trans_data[k1, c1_idx], trans_data[k2, c1_idx], markers.MarkerStyle.filled_markers[0], color = cmap(norm(0)), label=','.join(ct1) + ' Tumor', ms = 10, alpha=0.75)
                if c2_idx.shape[0] > 0:
                    ax.plot(trans_data[k1, c2_idx], trans_data[k2, c2_idx], markers.MarkerStyle.filled_markers[1], color = cmap(norm(2)), label=','.join(ct1) + ' Normal', ms = 10, alpha=0.75)
                if c3_idx.shape[0] > 0:
                    ax.plot(trans_data[k1, c3_idx], trans_data[k2, c3_idx], markers.MarkerStyle.filled_markers[2], color = cmap(norm(4)), label=','.join(ct2), ms = 10, alpha=0.75)

                ax.set_title('t-SNE Visualization for %s Events' % event_label_dict[event_type])
                ax.set_xticks(ax.get_xticks()[::2])
                ax.set_yticks(ax.get_yticks()[::2])
                ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
                ax.tick_params(axis = 'both', which = 'minor', labelsize = 10)
            #if k1 == (k - 1):
            #    ax.legend(numpoints = 1, ncol = 2, loc = 'center left', bbox_to_anchor = (1.2, 0.5))
            ax.legend(numpoints = 1, ncol = 1, loc = 'center left', bbox_to_anchor = (1.1, 0.6), frameon=False)
        plt.tight_layout()

        plt.savefig(os.path.join(plotdir, 'tsne_skl_%s.%s.%i.conf_%.2f.perp_%i%s.png' % (event_type, s, k, 1.0 - conf, perplexity, filter_tag)), dpi = 300, format = 'png', bbox_inches='tight')
        plt.savefig(os.path.join(plotdir, 'tsne_skl_%s.%s.%i.conf_%.2f.perp_%i%s.pdf' % (event_type, s, k, 1.0 - conf, perplexity, filter_tag)), dpi = 300, format = 'pdf', bbox_inches='tight')
        plt.close(fig)
    print('... done.')
 

