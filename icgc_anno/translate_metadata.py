import sys
import os
import cPickle
import gzip
import scipy as sp

META_REL = '/cluster/work/grlab/projects/ICGC/releaseTables/release_may2016.v1.4.tsv'
META_RNA = '/cluster/work/grlab/projects/ICGC/orig_data/metadata/per_aliquot_v2/rnaseq_metadata.histo.tsv'

def get_keys():
    keys = set()
    for line in open(META_REL, 'r'):
        for k in line.strip().split('\t'):
            keys.add(k)
        break
    for line in open(META_RNA, 'r'):
        for k in line.strip().split('\t'):
            keys.add(k)
        break
    return keys

def _handle_multi_entries(header, data):

    cols_of_interest = [['tumor_wgs_submitter_specimen_id', 
                         'tumor_wgs_icgc_specimen_id', 
                         'tumor_wgs_submitter_sample_id',
                         'tumor_wgs_icgc_sample_id',
                         'tumor_wgs_aliquot_id',
                         'tumor_wgs_oxog_score',
                         'tumor_wgs_ContEST',
                         'tumor_wgs_Stars',
                         'tumor_wgs_bwa_alignment_gnos_repo',
                         'tumor_wgs_bwa_alignment_gnos_id',
                         'is_mar2016_tumor_wgs_bwa_alignment',
                         'tumor_wgs_bwa_alignment_bam_file_name',
                         'tumor_wgs_minibam_gnos_repo',
                         'tumor_wgs_minibam_gnos_id',
                         'is_mar2016_tumor_wgs_minibam',
                         'tumor_wgs_minibam_bam_file_name',
                         'sanger_variant_calling_file_name_prefix',
                         'dkfz_embl_variant_calling_file_name_prefix',
                         'broad_variant_calling_file_name_prefix',
                         'muse_variant_calling_file_name_prefix',
                         'broad_tar_variant_calling_file_name_prefix',
                         'tumor_wgs_has_matched_rna_seq'],
                        ['tumor_rna_seq_submitter_specimen_id',
                         'tumor_rna_seq_icgc_specimen_id',
                         'tumor_rna_seq_submitter_sample_id',
                         'tumor_rna_seq_icgc_sample_id',
                         'tumor_rna_seq_aliquot_id',
                         'tumor_rna_seq_star_alignment_gnos_repo',
                         'tumor_rna_seq_star_alignment_gnos_id',
                         'is_mar2016_tumor_rna_seq_star_alignment',
                         'tumor_rna_seq_star_alignment_bam_file_name',
                         'tumor_rna_seq_tophat_alignment_gnos_repo',
                         'tumor_rna_seq_tophat_alignment_gnos_id',
                         'is_mar2016_tumor_rna_seq_tophat_alignment',
                         'tumor_rna_seq_tophat_alignment_bam_file_name']]

    for cols in cols_of_interest:
        c_idx = sp.where(sp.in1d(header, cols))[0]
        r_idx = sp.where([',' in x for x in data[:, c_idx[0]]])[0]

        for r in r_idx:
            data_ = sp.array([x.split(',') for x in data[r, c_idx]])
            assert len(data_.shape) > 1
            assert data_.shape[1] > 1

            for r2 in range(1, data_.shape[1]):
                data = sp.r_[data, data[r, :][sp.newaxis, :]]
                data[-1, c_idx] = data_[:, r2]
            data[r, c_idx] = data_[:, 0]

    return data
            


def _parse_metatable(fname):

    data = []
    for l, line in enumerate(open(fname, 'r')):
        sl = line.strip('\t').strip('\n').split('\t')
        if l == 0:
            header = sl
        else:
            data.append(sl)
    header = sp.array(header)
    data = sp.array(data)

    return (header, data)


def _merge_tables(header_rel, header_rna, data_rel, data_rna):

    aliquots_rel_t = data_rel[:, sp.where(header_rel == 'tumor_rna_seq_aliquot_id')[0][0]]
    aliquots_rel_n = data_rel[:, sp.where(header_rel == 'normal_rna_seq_aliquot_id')[0][0]]
    n_idx = sp.where(aliquots_rel_n != '')[0]
    aliquots_rel = sp.r_[aliquots_rel_t, aliquots_rel_n[n_idx]]
    data_rel = sp.r_[data_rel, data_rel[n_idx, :]]

    #k_idx = sp.where([',' in x for x in aliquots_rel])[0]
    #for k in k_idx:
    #    tmp = aliquots_rel[k].split(',')
    #    for j in range(1, len(tmp)):
    #        aliquots_rel = sp.append(aliquots_rel, tmp[j])
    #        data_rel = sp.r_[data_rel, data_rel[k, :][sp.newaxis, :]]
    #        data_rel[-1, sp.where(header_rel == 'tumor_rna_seq_aliquot_id')[0][0]] = tmp[j]
    #    data_rel[k, sp.where(header_rel == 'tumor_rna_seq_aliquot_id')[0][0]] = tmp[0]
    #    aliquots_rel[k] = tmp[0]

    aliquots_rna = data_rna[:, sp.where(header_rna == 'aliquot_id')[0][0]]
    a, b = sp.where(aliquots_rel == aliquots_rna[:, sp.newaxis])
    data_rna_ = sp.zeros((data_rel.shape[0], data_rna.shape[1]), dtype=data_rel.dtype)
    data_rna_[b, :] = data_rna[a, :]
    assert sp.all(data_rna_[b, sp.where(header_rna == 'icgc_donor_id')[0][0]] == data_rel[b, sp.where(header_rel == 'icgc_donor_id')[0][0]])
    
    k_idx = sp.where(~sp.in1d(header_rna, ['icgc_donor_id', 'wgs_exclusion_white_gray']))[0]
    header_rna = header_rna[k_idx]
    data_rna_ = data_rna_[:, k_idx]

    header = sp.r_[header_rel, header_rna]
    data = sp.c_[data_rel, data_rna_]
    del header_rel, header_rna, data_rel, data_rna, data_rna_

    return (header, data)


def translate(tr_pairs, first_only=False):
    '''Takes a list of pairs with keys that should be translated into each other

    Returns a list of dictionaries providing the requested translation
    '''

    metapickle = 'metadata.pickle.gz'
    if not os.path.exists(metapickle):
        (header_rel, data_rel) = _parse_metatable(META_REL)
        (header_rna, data_rna) = _parse_metatable(META_RNA)
        data_rel = _handle_multi_entries(header_rel, data_rel)

        ### merge the two tables into a single one
        (header, data) = _merge_tables(header_rel, header_rna, data_rel, data_rna)

        cPickle.dump((header, data), gzip.open(metapickle, 'w'), -1) 
    else:
        (header, data) = cPickle.load(gzip.open(metapickle, 'r'))

    dicts = []
    for source,target in tr_pairs:
        curr_dict = dict()
        idx1 = sp.where(header == source)[0][0]
        idx2 = sp.where(header == target)[0][0]
        s_idx = sp.argsort(data[:, idx1])
        data = data[s_idx, :]
        _, cnt = sp.unique(data[:, idx1], return_counts=True)
        cum = 0
        for c in cnt:
            if c == 1 or first_only:
                curr_dict[data[cum, idx1]] = data[cum, idx2]
            else:
                if sp.unique(data[cum:(cum+c), idx2]).shape[0] == 1:
                    curr_dict[data[cum, idx1]] = data[cum, idx2]
                else:
                    curr_dict[data[cum, idx1]] = ','.join(sp.unique(data[cum:(cum+c), idx2]))
            cum += c

        dicts.append(curr_dict)

    return dicts

    

if __name__ == "__main__":
    #translate([('normal_wgs_aliquot_id', 'icgc_donor_id')])
    print get_keys()
