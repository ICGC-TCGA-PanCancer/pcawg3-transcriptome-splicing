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

from paths import BASEDIR,BASEDIR_AS

event_types = ['mult_exon_skip', 'exon_skip'] #, 'intron_retention', 'alt_3prime', 'alt_5prime']
basedir_icgc = os.path.join(BASEDIR_AS, 'alternative_splicing')
metadata_icgc = os.path.join(BASEDIR, 'orig_data/metadata/per_aliquot_v2/rnaseq_metadata.tsv')
datadir = os.path.join(basedir_icgc, 'pca')
if not os.path.exists(datadir):
    os.makedirs(datadir)

if len(sys.argv) < 2:
    print('Usage: %s <k> [<conf> (0.01)]' % sys.argv[0], file=sys.stderr)
    sys.exit(1)
else:
    k = int(sys.argv[1])

if len(sys.argv) < 3:
    conf = 0.01
else:
    conf = float(sys.argv[2])

# create ICGC cancer type dictionary
(ct_dict, is_tumor_dict) = utils.get_ct_dict_metatable(metadata_icgc)

# create GTex type dictionary
sample_dict = {'BLDR':os.path.join(BASEDIR, 'gtex_tables/GTex_Bladder_samples.txt'),
               'BRST':os.path.join(BASEDIR, 'gtex_tables/GTex_Breast_samples.txt'),
               'KIDN':os.path.join(BASEDIR, 'gtex_tables/GTex_Kidney_samples.txt'),
               'LIVR':os.path.join(BASEDIR, 'gtex_tables/GTex_Liver_samples.txt'),
               'LUNG':os.path.join(BASEDIR, 'gtex_tables/GTex_Lung_samples.txt'),
               'PRST':os.path.join(BASEDIR, 'gtex_tables/GTex_Prostate_samples.txt'),
               'STOM':os.path.join(BASEDIR, 'gtex_tables/GTex_Stomach_samples.txt'),
               'THYR':os.path.join(BASEDIR, 'gtex_tables/GTex_Thyroid_samples.txt'),
               'UTER':os.path.join(BASEDIR, 'gtex_tables/GTex_Uterus_samples.txt'),
               'REM.':os.path.join(BASEDIR, 'gtex_tables/GTex_rest_wo_cells_samples.txt')}
gt_dict = utils.get_gt_dict(sample_dict)

for event_type in event_types:

    picklefile = '%s/pca2_%s.TN.conf_%.2f.pickle' % (datadir, event_type, 1.0 - conf)
    if True: #not os.path.exists(picklefile):
        ### load TCGA data from hdf5
        print('Loading data from TCGA hdf5')
        IN = h5py.File('%s/merge_graphs_%s_C3.counts.hdf5' % (basedir_icgc, event_type), 'r')
        c_idx = IN['conf_idx'][:].astype('int')
        #if event_type == 'exon_skip':
        #    a = IN['event_counts'][:, 4, :] + IN['event_counts'][:, 5, :]
        #    b = 2 * IN['event_counts'][:, 6, :]
        #elif event_type == 'intron_retention':
        #    a = IN['event_counts'][:, 1, :] # intron cov
        #    b = IN['event_counts'][:, 4, :] # intron conf
        #elif event_type in ['alt_3prime', 'alt_5prime']:
        #    a = IN['event_counts'][:, 3, :] # intron1 conf
        #    b = IN['event_counts'][:, 4, :] # intron2 conf
        strains = sp.array([x.split('.')[0] for x in IN['strains'][:]], dtype='str')

        ### get psi values
        psi = IN['psi'][:, c_idx]
        IN.close()

        ### only keep strains that are present in the metadata
        k_idx = sp.where(sp.in1d(strains, list(ct_dict.keys())))[0]
        strains = strains[k_idx]
        psi = psi[k_idx, :]

        ### sort by strains
        s_idx = sp.argsort(strains)
        strains = strains[s_idx]
        psi = psi[s_idx, :]

        ### get sample labels
        ctypes = sp.array([ct_dict[x] for x in strains])
        tn_labels = sp.array(['t' if is_tumor_dict[x] else 'n' for x in strains])

        ### subset to samples where we have normals
        ct_idx = sp.where(sp.in1d(ctypes, ctypes[sp.where(tn_labels == 'n')[0]]))[0]
        strains = strains[ct_idx]
        tn_labels = tn_labels[ct_idx]
        ctypes = ctypes[ct_idx]
        psi = psi[ct_idx, :]

        psi = psi.T

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
        #psi -= sp.mean(psi, axis=0)

        ### compute kernel
        print('Computing covariances')
        K = sp.cov([psi[i, :] for i in range(psi.shape[0])])

        ### PCA 
        print('Compute PCA ...')
        w_g, Vt_g = eigh(K)
        V_g = Vt_g.T
        w_g = w_g[::-1]
        V_g = V_g[::-1, :]
        print('... done')

        pickle.dump((w_g, V_g, ctypes, tn_labels, psi), open(picklefile, 'w'), -1)
    else:
        print('Loading data from pickle: %s' % picklefile)
        (w_g, V_g, ctypes, tn_labels, psi) = pickle.load(open(picklefile, 'r'))

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
            trans_data = sp.vstack([V_g[k1, :], V_g[k2, :]]).dot(psi)
            for idx, ct in enumerate(sp.unique(ctypes)):
                c_idx = sp.where(ctypes == ct)[0]
                if c_idx.shape[0] > 0:
                    ax.plot(trans_data[0, c_idx], trans_data[1, c_idx], markers.MarkerStyle.filled_markers[idx % 13], color = cmap(norm(idx)), label = ct, ms = 4, alpha=0.75)
            ax.set_title('PC %i vs %i' % (k1 + 1, k2 + 1))
            ax.set_xticks(ax.get_xticks()[::2])
            ax.set_yticks(ax.get_yticks()[::2])
            ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
            ax.tick_params(axis = 'both', which = 'minor', labelsize = 10)
        if k1 == (k - 1):
            ax.legend(numpoints = 1, ncol = 2, loc = 'center left', bbox_to_anchor = (1.2, 0.5))

    plt.tight_layout()
    plt.savefig('pca_%s.CT.%i.conf_%.2f.pdf' % (event_type, k, 1.0 - conf), dpi = 1200, format = 'pdf')
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
            trans_data = sp.dot(sp.vstack([V_g[k1, :], V_g[k2, :]]), psi)
            for idx, ct in enumerate(sp.unique(tn_labels)[::-1]):
                c_idx = sp.where(tn_labels == ct)[0]
                if c_idx.shape[0] > 0:
                    #ax.plot(V_g[k1, c_idx], V_g[k2, c_idx], 'o', color = cmap(norm(idx)), label = ct, ms = 4, alpha=0.75)
                    ax.plot(trans_data[0, c_idx], trans_data[1, c_idx], 'o', color = cmap(norm(idx)), label = ct, ms = 4, alpha=0.75)
            ax.set_title('PC %i vs %i' % (k1 + 1, k2 + 1))
            ax.set_xticks(ax.get_xticks()[::2])
            ax.set_yticks(ax.get_yticks()[::2])
            ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
            ax.tick_params(axis = 'both', which = 'minor', labelsize = 10)
        if k1 == (k - 1):
            ax.legend(numpoints = 1, ncol = 2, loc = 'center left', bbox_to_anchor = (1.2, 0.5))
    plt.tight_layout()

    plt.savefig('pca_%s.TN.%i.conf_%.2f.pdf' % (event_type, k, 1.0 - conf), dpi = 1200, format = 'pdf')
    print('... done.')

    sys.exit(0)

    print('Plotting individual ctype plots ...')
    sets = ['bladder', 'breast', 'lung', 'kidney']
    for s in sets:
        fig = plt.figure(figsize = (k*4, k*4))
        for k1 in range(0, k):
            cnt = 1
            for k2 in range(k1 + 1, k):
                ax = fig.add_subplot(k-1, k-1, (k1 * (k-1)) + cnt)
                cnt += 1
                trans_data = sp.vstack([V_g[k1, :], V_g[k2, :]]).dot(psi) 
                ### plot active elements
                if s == 'bladder':
                    norm = plt.Normalize(0, 6)
                    ct1 = 'BLCA'
                    ct2 = 'BLDR_GT'
                elif s == 'breast':
                    norm = plt.Normalize(0, 6)
                    ct1 = 'BRCA'
                    ct2 = 'BRST_GT'
                elif s == 'kidney':
                    norm = plt.Normalize(0, 6)
                    ct1 = 'KIRC'
                    ct2 = 'KIDN_GT'
                elif s == 'lung':
                    norm = plt.Normalize(0, 10)
                    ct1 = 'LUSC'
                    ct2 = 'LUNG_GT'
                    ct3 = 'LUAD'
                c1_idx = sp.where((ctypes == ct1) & (tn_labels == 't'))[0]
                c2_idx = sp.where((ctypes == ct1) & (tn_labels == 'n'))[0]
                c3_idx = sp.where(ctypes == ct2)[0]
                cc_idx = sp.hstack([c1_idx, c2_idx, c3_idx])
                if s == 'lung':
                    c4_idx = sp.where((ctypes == ct3) & (tn_labels == 't'))[0]
                    c5_idx = sp.where((ctypes == ct3) & (tn_labels == 'n'))[0]
                    cc_idx = sp.hstack([cc_idx, c4_idx, c5_idx])
                c6_idx = sp.where(~sp.in1d(sp.arange(ctypes.shape[0]), cc_idx))[0]
                ### plot inactive elements
                if c6_idx.shape[0] > 0:
                    ax.plot(trans_data[0, c6_idx], trans_data[1, c6_idx], markers.MarkerStyle.filled_markers[5], color = '0.6', label='other', ms = 4, alpha=0.5)
                ### plot active elements
                if c1_idx.shape[0] > 0:
                    ax.plot(trans_data[0, c1_idx], trans_data[1, c1_idx], markers.MarkerStyle.filled_markers[0], color = cmap(norm(0)), label=ct1 + '_T', ms = 4, alpha=0.75)
                if c2_idx.shape[0] > 0:
                    ax.plot(trans_data[0, c2_idx], trans_data[1, c2_idx], markers.MarkerStyle.filled_markers[1], color = cmap(norm(2)), label=ct1 + '_N', ms = 4, alpha=0.75)
                if c3_idx.shape[0] > 0:
                    ax.plot(trans_data[0, c3_idx], trans_data[1, c3_idx], markers.MarkerStyle.filled_markers[2], color = cmap(norm(4)), label=ct2, ms = 4, alpha=0.75)
                if s == 'lung':
                    if c4_idx.shape[0] > 0:
                        ax.plot(trans_data[0, c4_idx], trans_data[1, c4_idx], markers.MarkerStyle.filled_markers[3], color = cmap(norm(6)), label=ct3 + '_T', ms = 4, alpha=0.75)
                    if c4_idx.shape[0] > 0:
                        ax.plot(trans_data[0, c5_idx], trans_data[1, c5_idx], markers.MarkerStyle.filled_markers[4], color = cmap(norm(8)), label=ct3 + '_N', ms = 4, alpha=0.75)

                ax.set_title('PC %i vs %i' % (k1 + 1, k2 + 1))
                ax.set_xticks(ax.get_xticks()[::2])
                ax.set_yticks(ax.get_yticks()[::2])
                ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
                ax.tick_params(axis = 'both', which = 'minor', labelsize = 10)
            if k1 == (k - 1):
                ax.legend(numpoints = 1, ncol = 2, loc = 'center left', bbox_to_anchor = (1.2, 0.5))
        plt.tight_layout()

        plt.savefig('pca_%s.%s.%i.conf_%.2f.pdf' % (event_type, s, k, 1.0 - conf), dpi = 1200, format = 'pdf')
    print('... done.')
 

