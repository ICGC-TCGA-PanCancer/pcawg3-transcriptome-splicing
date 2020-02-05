import scipy as sp
import scipy.io as scio
import sys
import h5py
import os
import pickle
import glob
import pdb
import re
sys.path.append('/cluster/home/akahles/git/software/spladder/python')
import modules.classes.gene

from samples import intron_sample, psi_sample, event_sample

def _load_tsv(fname):
    
    data = []
    header = []
    for i, line in enumerate(open(fname, 'r')):
        sl = line.strip().split('\t')
        if i == 0:
            header = sl
        else:
            data.append(sl)
    data = sp.array(data, dtype='str')
    header = sp.array(header, dtype='str')

    return (data, header)


def get_ct_dict(sample_pattern):

    ### construct cancer type dictionary
    ct_data = sp.array([x.split('/')[-2:] for x in glob.glob(sample_pattern)])
    ct_data[:, 1] = sp.array(['-'.join(x.split('-')[:4]) for x in ct_data[:, 1]])

    return  dict((ct_data[i, 1], ct_data[i, 0]) for i in range(ct_data.shape[0]))


def get_ct_dict_metatable(fname, style='icgc'):

    ### load metadata table filename
    (data, header) = _load_tsv(fname)

    if style == 'icgc':
        a_idx = sp.where(header == 'analysis_id')[0][0]
        c_idx = sp.where(header == 'project_code')[0][0]
        is_tumor_idx = sp.where(header == 'is_tumour')[0][0]

        ct_dict = dict((data[i, a_idx], data[i, c_idx]) for i in range(data.shape[0])) 
        is_tumor_dict = dict((data[i, a_idx], data[i, is_tumor_idx] == 'yes') for i in range(data.shape[0]))
    elif style == 'pancan_rerun':
        a_idx = sp.where(header == 'analysis_id')[0][0]
        c_idx = sp.where(header == 'study')[0][0]
        is_tumor_idx = sp.where(header == 'is_normal')[0][0]
        
        ct_dict = dict((data[i, a_idx], data[i, c_idx]) for i in range(data.shape[0])) 
        is_tumor_dict = dict((data[i, a_idx], data[i, is_tumor_idx] == 'False') for i in range(data.shape[0]))
    
    return (ct_dict, is_tumor_dict)

def get_ht_dict_metatable(fname_metadata, fname_histotypes):
    
    ### load metadata table
    (metadata, metaheader) = _load_tsv(fname_metadata)

    ### load histotype data table
    (histodata, histoheader) = _load_tsv(fname_histotypes)

    ### specimen id in histotable
    sh_idx = sp.where(histoheader == '# icgc_specimen_id')[0][0]
    ### specimen id in metatable
    sm_idx = sp.where(metaheader == 'icgc_specimen_id')[0][0]
    ### histotype id
    h_idx = sp.where(histoheader == 'histology_abbreviation')[0][0]
    ### analysis id
    a_idx = sp.where(metaheader == 'analysis_id')[0][0]

    histo = dict()
    for i in range(histodata.shape[0]):
        sh = histodata[i, sh_idx]
        h = histodata[i, h_idx]
        if sh in histo:
            assert histo[sh] == h
        else:
            histo[sh] = h

    aid2histo = dict()
    for i in range(metadata.shape[0]):
        sm = metadata[i, sm_idx]
        aid = metadata[i, a_idx]
        assert not aid in aid2histo
        try:
            aid2histo[aid] = histo[sm]
        except KeyError:
            aid2histo[aid] = 'NA'

    return aid2histo


def get_gt_dict(sample_dict):

    ### construct gtex tissue type dictionary
    gt_dict = dict()
    for tt in sample_dict:
        for line in open(sample_dict[tt], 'r'):
            gt_dict[line.strip()] = tt

    return gt_dict
            

def get_expression_files(fnames, rules=None):

    ### load gene expression counts from hdf5
    for f, fname in enumerate(fnames):
        IN = h5py.File(fname, 'r')
        if f == 0:
            exp_cnt = IN['counts'][:]
            exp_gids = IN['gids'][:]
            exp_sids = IN['sids'][:]
            if rules is not None:
                is_tumor = sp.array([rules[f].match(x) is not None for x in exp_sids], dtype='bool') 
        else:
            exp_cnt = sp.c_[exp_cnt, IN['counts'][:]]
            exp_sids = sp.r_[exp_sids, IN['sids'][:]]
            if rules is not None:
                is_tumor = sp.r_[is_tumor, sp.array([rules[f].match(x) is not None for x in IN['sids'][:]], dtype='bool')]
            assert(exp_gids.shape[0] == IN['gids'].shape[0])
        IN.close()

    if rules is not None:
        return (exp_cnt, exp_gids, exp_sids, is_tumor)
    else:
        return (exp_cnt, exp_gids, exp_sids)


def get_iso_counts(file, event_type, idx=None, normalize=True):
    """Loads isoform data for event_type from event hdf5 file"""

    print("loading %s ..." % file, file=sys.stderr)
    IN = h5py.File(file, 'r')
    if idx is not None:
        e_idx = idx
    else:
        e_idx = sp.arange(IN['event_counts'].shape[0])

    ### get normalization factor
    if normalize:
        if event_type in ['exon_skip', 'intron_retention']:
            sf = sp.sum(IN['event_counts'][:, 2, :], axis=0) + sp.sum(IN['event_counts'][:, 3, :], axis=0)
        elif event_type in ['alt_3prime', 'alt_5prime']:
            sf = sp.sum(IN['event_counts'][:, 2, :], axis=0)
        sf /= sp.median(sf)
    else:
        sf = 1.0

    if event_type == 'exon_skip':
        ### exon is included
        iso1 = sp.minimum(IN['event_counts'][:, 4, :], IN['event_counts'][:, 5, :]).astype('float') * sf
        ### exon is skipped
        iso2 = IN['event_counts'][:, 6, :].astype('float') * sf
    elif event_type == 'intron_retention':
        ### intron is retained
        iso1 = IN['event_counts'][:, 1, :].astype('float') * sf
        ### intron is spliced
        iso2 = IN['event_counts'][:, 4, :].astype('float') * sf
    elif event_type in ['alt_3prime', 'alt_5prime']:
        ### intron1 is used
        iso1 = IN['event_counts'][:, 3, :].astype('float') * sf
        ### intron2 is used
        iso2 = IN['event_counts'][:, 4, :].astype('float') * sf
    eid = sp.array(['%s.' % IN['gene_chr'][IN['gene_idx'][i].astype('int')] + '.'.join(IN['event_pos'][i, :].astype('str')) for i in range(IN['gene_idx'].shape[0])])
    gid = sp.array([IN['gene_names'][IN['gene_idx'][i].astype('int')] for i in range(IN['gene_idx'].shape[0])])
    strains = IN['strains'][:]
    IN.close()

    return (iso1, iso2, eid, gid, strains)

def get_psi_counts(file, event_type, idx=None, annotation=None):
    """Loads psi data for event_type from event hdf5 file"""

    print("loading %s ..." % file, file=sys.stderr)
    IN = h5py.File(file, 'r')
    strains = IN['strains'][:].astype('str')
    
    ### if event_type is defined, we continue here 
    if idx is not None:
        e_idx = idx
    else:
        e_idx = sp.arange(IN['event_counts'].shape[0])

    ### extract introns from single events
    psi = IN['psi'][:].T
    if event_type == 'mutex_exons':
        eid = sp.array(['%s.' % IN['gene_chr'][IN['gene_idx'][i].astype('int')] + '.'.join(IN['event_pos'][i, :].ravel().astype('str')) for i in range(IN['gene_idx'].shape[0])])
    else:
        eid = sp.array(['%s.' % IN['gene_chr'][IN['gene_idx'][i].astype('int')] + '.'.join(IN['event_pos'][i, :].astype('str')) for i in range(IN['gene_idx'].shape[0])])
       # tmp1 = sp.array(['%s.' % IN['gene_chr'][IN['gene_idx'][i].astype('int')] + '.'.join(IN['event_pos'][i, [1, 2]].astype('str')) for i in range(IN['gene_idx'].shape[0])])
       # tmp2 = sp.array(['%s.' % IN['gene_chr'][IN['gene_idx'][i].astype('int')] + '.'.join(IN['event_pos'][i, [3, 4]].astype('str')) for i in range(IN['gene_idx'].shape[0])])
       # tmp3 = sp.array(['%s.' % IN['gene_chr'][IN['gene_idx'][i].astype('int')] + '.'.join(IN['event_pos'][i, [1, 4]].astype('str')) for i in range(IN['gene_idx'].shape[0])])
    gid = sp.array([IN['gene_names'][IN['gene_idx'][i].astype('int')] for i in range(IN['gene_idx'].shape[0])])
       # eid = sp.c_[tmp1, tmp2, tmp3].T
       # gid = sp.c_[gid_, gid_, gid_].T
   # elif event_type in ['alt_3prime', 'alt_5prime']:
   #     coords = sp.c_[IN['event_pos'][:, [1, 2]], IN['event_pos'][:, [1,4]], IN['event_pos'][:, [3, 4]]].reshape(IN['event_pos'].shape[1] * 3, 2)
   #     tmp1 = sp.array(['%s.' % IN['gene_chr'][IN['gene_idx'][i].astype('int')] + '.'.join(IN['event_pos'][i, [1, 2]].astype('str')) for i in range(IN['gene_idx'].shape[0])])
   #     tmp2 = sp.array(['%s.' % IN['gene_chr'][IN['gene_idx'][i].astype('int')] + '.'.join(IN['event_pos'][i, [1, 4]].astype('str')) for i in range(IN['gene_idx'].shape[0])])
   #     tmp3 = sp.array(['%s.' % IN['gene_chr'][IN['gene_idx'][i].astype('int')] + '.'.join(IN['event_pos'][i, [3, 4]].astype('str')) for i in range(IN['gene_idx'].shape[0])])
   #     gid_ = sp.array([IN['gene_names'][IN['gene_idx'][i].astype('int')] for i in range(IN['gene_idx'].shape[0])])
   #     gid = sp.c_[gid_, gid_, gid_].reshape(gid_.shape[0] * 3,)
   #     eid = sp.c_[tmp1, tmp2, tmp3].reshape(tmp1.shape[0] * 3,)
   #     k_idx = sp.where(coords[:, 1] - coords[:, 0] > 0)[0]
   #     eid = eid[k_idx]
   #     gid = gid[k_idx]
   #     assert(eid.shape[0] == gid.shape[0])
   # elif event_type == 'intron_retention':
   #     eid = sp.array(['%s.' % IN['gene_chr'][IN['gene_idx'][i].astype('int')] + '.'.join(IN['event_pos'][i, [1, 2]].astype('str')) for i in range(IN['gene_idx'].shape[0])])
   #     gid = sp.array([IN['gene_names'][IN['gene_idx'][i].astype('int')] for i in range(IN['gene_idx'].shape[0])])
    IN.close()

    ### make events unique
    s_idx = sp.argsort(eid)
    eid = eid[s_idx]
    gid = gid[s_idx]
    psi = psi[s_idx, :]
    eid, u_idx = sp.unique(eid, return_index=True)
    gid = gid[u_idx]
    psi = psi[u_idx, :]

    return (psi, eid, gid, strains)


def get_intron_counts(file, event_type, idx=None, normalize=True, annotation=None):
    """Loads intron data for event_type from event hdf5 file"""

    print("loading %s ..." % file, file=sys.stderr)
    IN = h5py.File(file, 'r')
    strains = IN['strains'][:].astype('str')
    
    ### if event_type is None take all introns that we ever counted
    if event_type is None:
        genes = pickle.load(open(annotation, 'r'))[0]
        
        introns = IN['edges'][:]
        gids = IN['gene_names'][:]
        gid = []
        eid = []
        ### iterate over edges we have found in in the counts
        for i,j in enumerate(IN['gene_ids_edges'][:, 0]):
            if i > 0 and i % 10 == 0:
                print('.', end=' ')
                if i % 100 == 0:
                    print('%i / %i' % (i, IN['gene_ids_edges'].shape[0]))

            g_idx = int(IN['gene_ids_edges'][i, 0])
            gid.append(gids[g_idx])
            assert(genes[g_idx].name == gid[-1])
            try:
                start, end = sp.unravel_index(int(IN['edge_idx'][i]), genes[g_idx].segmentgraph.seg_edges.shape)
            except:
                import pdb
                pdb.set_trace()
            eid.append('%s.%i.%i' % (genes[g_idx].chr, genes[g_idx].segmentgraph.segments[1, start], genes[g_idx].segmentgraph.segments[0, end]))

        gid = sp.array(gid, dtype='str')
        eid = sp.array(eid, dtype='str')

        if strains[0].startswith('hg19'):
            strains = sp.array([x.split(':')[1].split('.')[0] for x in strains], dtype='str')
        else:
            strains = sp.array(strains, dtype='str')
        return (introns, eid, gid, strains)
                
    ### if event_type is defined, we continue here 
    if idx is not None:
        e_idx = idx
    else:
        e_idx = sp.arange(IN['event_counts'].shape[0])

    ### get normalization factor
    if normalize:
        if event_type in ['exon_skip', 'intron_retention']:
            sf = sp.sum(IN['event_counts'][:, 2, :], axis=0) + sp.sum(IN['event_counts'][:, 3, :], axis=0)
        elif event_type in ['alt_3prime', 'alt_5prime']:
            sf = sp.sum(IN['event_counts'][:, 2, :], axis=0)
        sf /= sp.median(sf)
    else:
        sf = 1.0

    ### extract introns from single events
    if event_type == 'exon_skip':
        introns = sp.r_[IN['event_counts'][:, 4, :], IN['event_counts'][:, 5, :], IN['event_counts'][:, 6, :]].astype('float') * sf
        tmp1 = sp.array(['%s.' % IN['gene_chr'][IN['gene_idx'][i].astype('int')] + '.'.join(IN['event_pos'][1:3, i].astype('str')) for i in range(IN['gene_idx'].shape[0])])
        tmp2 = sp.array(['%s.' % IN['gene_chr'][IN['gene_idx'][i].astype('int')] + '.'.join(IN['event_pos'][3:5, i].astype('str')) for i in range(IN['gene_idx'].shape[0])])
        tmp3 = sp.array(['%s.' % IN['gene_chr'][IN['gene_idx'][i].astype('int')] + '.'.join(IN['event_pos'][[1,4], i].astype('str')) for i in range(IN['gene_idx'].shape[0])])
        gid_ = sp.array([IN['gene_names'][IN['gene_idx'][i].astype('int')] for i in range(IN['gene_idx'].shape[0])])
        eid = sp.r_[tmp1, tmp2, tmp3]
        gid = sp.r_[gid_, gid_, gid_]
    elif event_type in ['alt_3prime', 'alt_5prime']:
        introns = sp.zeros((IN['event_counts'].shape[0] * 2, IN['event_counts'].shape[2]), dtype='float')
        introns[::2, :] = IN['event_counts'][:, 3, :].astype('float') * sf
        introns[1::2, :] = IN['event_counts'][:, 4, :].astype('float') * sf
        coords = sp.c_[IN['event_pos'][1:3, :].T, IN['event_pos'][[1,4], :].T, IN['event_pos'][3:5, :].T].reshape(IN['event_pos'].shape[1] * 3, 2)
        tmp1 = sp.array(['%s.' % IN['gene_chr'][IN['gene_idx'][i].astype('int')] + '.'.join(IN['event_pos'][1:3, i].astype('str')) for i in range(IN['gene_idx'].shape[0])])
        tmp2 = sp.array(['%s.' % IN['gene_chr'][IN['gene_idx'][i].astype('int')] + '.'.join(IN['event_pos'][[1,4], i].astype('str')) for i in range(IN['gene_idx'].shape[0])])
        tmp3 = sp.array(['%s.' % IN['gene_chr'][IN['gene_idx'][i].astype('int')] + '.'.join(IN['event_pos'][3:5, i].astype('str')) for i in range(IN['gene_idx'].shape[0])])
        gid_ = sp.array([IN['gene_names'][IN['gene_idx'][i].astype('int')] for i in range(IN['gene_idx'].shape[0])])
        gid = sp.c_[gid_, gid_, gid_].reshape(gid_.shape[0] * 3,)
        eid = sp.c_[tmp1, tmp2, tmp3].reshape(tmp1.shape[0] * 3,)
        k_idx = sp.where(coords[:, 1] - coords[:, 0] > 0)[0]
        eid = eid[k_idx]
        gid = gid[k_idx]
        assert(eid.shape[0] == introns.shape[0])
        assert(eid.shape[0] == gid.shape[0])
    elif event_type == 'intron_retention':
        introns = IN['event_counts'][:, 4, :].astype('float') * sf
        eid = sp.array(['%s.' % IN['gene_chr'][IN['gene_idx'][i].astype('int')] + '.'.join(IN['event_pos'][1:3, i].astype('str')) for i in range(IN['gene_idx'].shape[0])])
        gid = sp.array([IN['gene_names'][IN['gene_idx'][i].astype('int')] for i in range(IN['gene_idx'].shape[0])])

    IN.close()

    ### make introns unique
    s_idx = sp.argsort(eid)
    eid = eid[s_idx]
    gid = gid[s_idx]
    introns = introns[s_idx, :]
    eid, u_idx = sp.unique(eid, return_index=True)
    gid = gid[u_idx]
    introns = introns[u_idx, :]

    return (introns, eid, gid, strains)


def preprocess_counts(event_type, basedir, hdf5_files, rules=None):

    samples = dict()

    for key in list(hdf5_files.keys()):
        if event_type is None:
            picklefile = os.path.join(basedir, 'comparisons', 'pickle', 'compare_TCGAvsNormals.all_events.%s.hdf5' % (key))
        else:
            picklefile = os.path.join(basedir, 'comparisons', 'pickle', 'compare_TCGAvsNormals.%s.%s.hdf5' % (event_type, key))
        if not os.path.exists(os.path.join(basedir, 'comparisons', 'pickle')):
            os.makedirs(os.path.join(basedir, 'comparisons', 'pickle'))

        ### load event counts from count hdf5 if processed hdf5 is not available
        if not os.path.exists(picklefile):

            (_iso1, _iso2, _eid, _gid, _strains) = get_iso_counts(hdf5_files[key], event_type, normalize=False)

            samples[key] = event_sample(_iso1.T, _iso2.T, _eid, _gid, _strains)
            samples[key].to_hdf5(picklefile)
        ### load event counts from preprocessed hdf5
        else:
            samples[key] = event_sample.from_hdf5(picklefile)

        samples[key].strains = sp.array([re.sub(r'.aligned$', '', x) for x in samples[key].strains], dtype='str') 

        ### subset samples according to rules if present
        if rules is not None and key in rules:
            samples[key].subset_by_strain(rules[key])

    ### assert that the files we use quantified the same events ...
    for i,k in enumerate(samples):
        if i == 0:
            lastk = k
            continue
        assert(sp.all(samples[k].eid == samples[lastk].eid))
        assert(sp.all(samples[k].gid == samples[lastk].gid))
        lastk = k


    return samples

def preprocess_psi_counts(event_type, basedir, hdf5_files, rules=None, annotation=None, project=None):

    samples = dict()

    for key in list(hdf5_files.keys()):
        if event_type is None:
            picklefile = os.path.join(basedir, 'comparisons', 'pickle', 'compare_TCGAvsNormals.all_psi.%s.hdf5' % (key))
        else:
            picklefile = os.path.join(basedir, 'comparisons', 'pickle', 'compare_TCGAvsNormals.psi.%s.%s.hdf5' % (event_type, key))
        if not os.path.exists(os.path.join(basedir, 'comparisons', 'pickle')):
            os.makedirs(os.path.join(basedir, 'comparisons', 'pickle'))

        ### load psi values from count hdf5 if psi hdf5 is not available
        if not os.path.exists(picklefile):

            (_psi, _eid, _gid, _strains) = get_psi_counts(hdf5_files[key], event_type, annotation=annotation)

            samples[key] = psi_sample(_psi, _eid, _gid, _strains)
            samples[key].to_hdf5(picklefile)

        ### load psi from available hdf5
        else:
            samples[key] = psi_sample.from_hdf5(picklefile)

        ### TODO: Make this more general
        if project is None:
            samples[key].strains = sp.array(['.'.join(x.split('.')[:1]) for x in samples[key].strains], dtype='str')
        elif project == 'pancan':
            samples[key].strains = sp.array(['.'.join(x.split('.')[:2]) for x in samples[key].strains], dtype='str')

        ### subset samples according to rules if present
        if rules is not None and key in rules:
            samples[key].subset_by_strain(rules[key])

    ### assert that the files we use quantified the same events ...
    for i,k in enumerate(samples):
        if i == 0:
            lastk = k
            continue
        assert(sp.all(samples[k].eid == samples[lastk].eid))
        assert(sp.all(samples[k].gid == samples[lastk].gid))
        lastk = k

    return samples


def preprocess_intron_counts(event_type, basedir, hdf5_files, rules=None, annotation=None, dataset='icgc', tag=''):

    samples = dict()

    for key in list(hdf5_files.keys()):
        if event_type is None:
            picklefile = os.path.join(basedir, 'comparisons', 'pickle', 'compare_TCGAvsNormals.all_introns%s.%s.hdf5' % (tag, key))
        else:
            picklefile = os.path.join(basedir, 'comparisons', 'pickle', 'compare_TCGAvsNormals.introns%s.%s.%s.hdf5' % (tag, event_type, key))
        if not os.path.exists(os.path.join(basedir, 'comparisons', 'pickle')):
            os.makedirs(os.path.join(basedir, 'comparisons', 'pickle'))

        ### load introns from hdf5 is pickle not available
        if not os.path.exists(picklefile):

            (_introns, _eid, _gid, _strains) = get_intron_counts(hdf5_files[key], event_type, normalize=False, annotation=annotation)

            samples[key] = intron_sample(_introns, _eid, _gid, _strains)
            samples[key].to_hdf5(picklefile)

        ### load introns from available hdf5
        else:
            samples[key] = intron_sample.from_hdf5(picklefile)

        ### TODO: Make this more general
        if dataset == 'pancanatlas':
            samples[key].strains = sp.array(['.'.join(x.split('.')[:2]) for x in samples[key].strains], dtype='str')
        else:
            samples[key].strains = sp.array(['.'.join(x.split('.')[:1]) for x in samples[key].strains], dtype='str')

        ### subset samples according to rules if present
        if rules is not None and key in rules:
            samples[key].subset_by_strain(rules[key])

    ### assert that the files we use quantified the same events ...
    for i,k in enumerate(samples):
        if i == 0:
            lastk = k
            continue
        assert(sp.all(samples[k].eid == samples[lastk].eid))
        assert(sp.all(samples[k].gid == samples[lastk].gid))
        lastk = k

    return samples


def gen_plot_directory(options):
    """Assemble output directory name for plotting.
    
    From the information in options the appropriate name of
    the plotting output directpry is built. If the directory 
    should not exist yet, it will be created.
    
    Returns:
        name of the plotting directory
    """

    ### add thresholds to directory name
    if options.thresholds is not None:
        plotdir = os.path.join(options.plotdir, '_'.join(['%s%i' % (x, y) for x, y in options.thresholds.items()]))
    else:
        plotdir = os.path.join(options.plotdir, 'no_thresh')

    ### add expression thresholds to directory name
    if options.max_norm_expression_frac is not None:
        if type(options.max_norm_expression_frac) is not dict:
            plotdir = '%s_ex-all%.2f' % (plotdir, options.max_norm_expression_frac)
        else:
            plotdir = '%s_%s' % (plotdir, '_'.join(['ex-%s%.2f' % (x, y) for x, y in options.max_norm_expression_frac.items()]))
    else:
        plotdir = '%s_ex-no' % (plotdir)

    ### add testing strategy to directory name
    plotdir = '%s.%s' % (plotdir, options.strategy)
        
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)

    return plotdir


def get_anno_from_mat(anno_fname):

    ### load annotation
    anno = scio.loadmat(anno_fname)
    anno = anno['genes'][0, :] 

    ### generate list of annotation IDs of introns
    anno_ids = []
    for gene in anno:
        chrm = gene['chr'][0]
        for t in range(gene['exons'].shape[1]):
            if gene['exons'][0, t].shape[0] > 1:
                anno_ids.extend(['%s.%i.0.%i.0' % (chrm, gene['exons'][0, t][i-1, 1], gene['exons'][0, t][i, 0]) for i in range(1, gene['exons'][0, t].shape[0])])
    anno_ids = sp.array([x.replace('chr', '') for x in anno_ids], dtype='str')

    return sp.unique(sp.sort(anno_ids))


