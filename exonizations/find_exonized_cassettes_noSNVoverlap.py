import sys
import os
import scipy as sp
import scipy.sparse as spsp
import numpy.random as npr
npr.seed(23)
import pickle
import re
import intervaltree as it
import h5py
import scipy.stats as spst

sys.path.append('/cluster/home/akahles/git/software/spladder/python')
import modules.utils as mu

sys.path.append('../icgc_anno')
import translate_metadata as tm

sys.path.append('../icgc_utils')
import names as un
lookup = un.get_lookup_complete()

sys.path.append('..')
from paths import BASEDIR,BASEDIR_AS

event_type = 'exon_skip'
event_in = os.path.join(BASEDIR_AS, 'alternative_splicing', 'merge_graphs_%s_C2.pickle' % event_type)
candidate_out_step1 = os.path.join(BASEDIR_AS, 'alternative_splicing', 'merge_graphs_%s_C2.exonize_candidates.pickle' % event_type)
candidate_out_step2 = os.path.join(BASEDIR_AS, 'alternative_splicing', 'merge_graphs_%s_C2.exonize_candidates_step2.pickle' % event_type)
candidate_out_step3 = os.path.join(BASEDIR_AS, 'alternative_splicing', 'merge_graphs_%s_C2.exonize_candidates_step3.pickle' % event_type)
candidate_out_step4 = os.path.join(BASEDIR_AS, 'alternative_splicing', 'merge_graphs_%s_C2.exonize_candidates_step4.pickle' % event_type)
recomputed_psi = os.path.join(BASEDIR_AS, 'alternative_splicing', 'merge_graphs_%s_C2.psi_recomputed.pickle' % event_type)

VARIANTS = os.path.join(BASEDIR, 'qtl_analysis/variants/mccalls/October_2016_whitelist_2583.snv_mnv_indel.sorted.sparse.hdf5')
coding_genes = sp.loadtxt(os.path.join(BASEDIR, 'annotation/gencode.v19.annotation.hs37d5_chr.gtf.coding_genes.txt'), delimiter='\t', dtype='str')     
read_thresh = 3

### prepare bam file dict
bam_dict = dict()
for line in open(os.path.join(BASEDIR_AS, 'alternative_splicing', 'sample_list_merged.txt'), 'r'):
    bam_id = re.sub(r'.bam$', '', line.strip().split('/')[-1])
    bam_dict[bam_id] = line.strip()

### get ID mappings
[strain_dict, tumor_dict, primary_dict, project_dict] = tm.translate([('analysis_id', 'icgc_donor_id'), ('analysis_id', 'is_tumour'), ('analysis_id', 'specimen_type'), ('icgc_donor_id', 'project_code')])

### get conf_idx
EV = h5py.File(os.path.join(BASEDIR_AS, 'alternative_splicing', 'merge_graphs_%s_C2.counts.hdf5' % event_type), 'r')

conf_idx = EV['conf_idx'][:]
strains_ev = sp.array([x.split('.')[0] for x in EV['strains'][:]])
files_ev = sp.array([bam_dict[x] for x in EV['strains'][:]]) 
donors_ev = sp.array([strain_dict[x] if x in strain_dict else 'NA' for x in strains_ev])

### only keep tumor samples
st_idx = sp.where([tumor_dict[x] == 'yes' if x in tumor_dict else False for x in strains_ev])[0]  
### only keep primary samples
k_idx = sp.where([primary_dict[x].startswith('Primary') if x in primary_dict else 'False' for x in strains_ev[st_idx]])[0]
st_idx = st_idx[k_idx]
_, u_idx = sp.unique(strains_ev[st_idx], return_index=True)
st_idx = st_idx[u_idx]

strains_ev = sp.array([strain_dict[x] if x in strain_dict else 'NA' for x in strains_ev])

### load event data
print('loading event data ...')
data = pickle.load(open(event_in, 'r'))

is_candidate = []
if not os.path.exists(candidate_out_step1):
    ### load gene annotations
    print('loading genes ...')
    genes, _ = pickle.load(open(os.path.join(BASEDIR_AS, 'alternative_splicing', 'spladder', 'genes_graph_conf2.merge_graphs.pickle'), 'rb'), encoding='latin')
    ### generate intron lists per gene
    introns = []
    anno_exons = dict()
    chrms = sp.array([_.chr for _ in genes])
    for chrm in sp.unique(chrms):
        anno_exons[chrm] = it.IntervalTree() 
    for i in range(genes.shape[0]):
        curr_introns = sp.zeros((0, 2), dtype='int')
        for j in range(len(genes[i].exons)):
            exons = genes[i].exons[j]
            if exons.shape[0] > 1:
                exons = mu.sort_rows(exons)
                curr_introns = sp.r_[curr_introns, exons.ravel()[1:-1].reshape(exons.shape[0] - 1, 2)]
            for k in range(exons.shape[0]):
                anno_exons[genes[i].chr][exons[k, 0]:exons[k, 1]] = i
        introns.append(mu.unique_rows(curr_introns))

    ### iterate over events and flag them:
    ###     - if they contain an intron that can not be found in the annotation
    ###     - if the inner exon overlaps an annotated exon
    cnt = 0
    for i in range(data.shape[0]):
        if i > 0 and i % 1000 == 0:
            sys.stdout.write('.')
            if i % 10000 == 0:
                sys.stdout.write('%i/%i - %i candidates\n' % (i, data.shape[0], cnt))
            sys.stdout.flush()
        ### check whether skipping is annotated
        if not mu.ismember([data[i].exons1[0, 1], data[i].exons1[1, 0]], introns[data[i].gene_idx], rows=True):
            is_candidate.append(0)
            continue
        ### check whether there exists an annotated exon overlapping to the inner exon
        tmp = 0
        for itv in anno_exons[data[i].chr].search(data[i].exons2[1, 0], data[i].exons2[1, 1]):
            tmp += 1
        if tmp > 0:
            is_candidate.append(0)
            continue
        ### check whether both flankings are not annotated
        if mu.ismember([data[i].exons2[0, 1], data[i].exons2[1, 0]], introns[data[i].gene_idx], rows=True) or \
           mu.ismember([data[i].exons2[1, 1], data[i].exons2[2, 0]], introns[data[i].gene_idx], rows=True):
            is_candidate.append(0)
            continue
        is_candidate.append(1)
        cnt += 1
    is_candidate = sp.array(is_candidate, dtype='bool')
    pickle.dump(is_candidate, open(candidate_out_step1, 'wb'), -1)
else:
    print('%s already exists' % candidate_out_step1)
    is_candidate = pickle.load(open(candidate_out_step1, 'rb'), encoding='latin')
print('considering %i events' % sp.sum(is_candidate))


### build GTEx intron index
### check which novel exons we would also find in GTEx
print('building GTEx intron index', file=sys.stderr)
GTEX = os.path.join(BASEDIR_AS, 'alternative_splicing_gtex/spladder/genes_graph_conf2.merge_graphs.pickle')
GTEX_INTRONS = os.path.join(BASEDIR_AS, 'alternative_splicing_gtex/spladder/genes_graph_conf2.merge_graphs.introns.pickle')
if not os.path.exists(GTEX_INTRONS):
    print('generating GTEx intron list')
    genes = pickle.load(open(GTEX, 'rb'), encoding='latin')
    gtex_introns = dict()
    gtex_exons = dict()

    for j,g in enumerate(genes[0]):
        if j > 0 and j % 100 == 0:
            sys.stdout.write('.')
            if j % 1000 == 0:
                sys.stdout.write('%i/%i\n' % (j, genes[0].shape[0]))
            sys.stdout.flush()
        sg_edges = spsp.coo_matrix((g.splicegraph_edges_data, (g.splicegraph_edges_row, g.splicegraph_edges_col)), g.splicegraph_edges_shape).toarray()
        x, y = sp.where(sp.triu(sg_edges))
        for i in range(x.shape[0]):
            intron = ':'.join([g.chr, str(g.splicegraph.vertices[1, x[i]]), str(g.splicegraph.vertices[0, y[i]])])
            try:
                gtex_introns[g.chr].append(intron)
            except KeyError:
                gtex_introns[g.chr] = [intron]
        for i in range(g.splicegraph.vertices.shape[1]):
            exon = ':'.join([g.chr, str(g.splicegraph.vertices[0, i]), str(g.splicegraph.vertices[1, i])])
            try:
                gtex_exons[g.chr].append(exon)
            except KeyError:
                gtex_exons[g.chr] = [exon]

    for chrm in gtex_introns:
        gtex_introns[chrm] = set(gtex_introns[chrm])
        gtex_exons[chrm] = set(gtex_exons[chrm])
    
    pickle.dump((gtex_introns, gtex_exons), open(GTEX_INTRONS, 'wb'), -1)
else:
    print('loading GTEx intron list from %s' % GTEX_INTRONS)
    (gtex_introns, gtex_exons) = pickle.load(open(GTEX_INTRONS, 'rb'), encoding='latin')

if not os.path.exists(candidate_out_step2):
    keep_idx = []
    for ii,i in enumerate(sp.where(is_candidate)[0]):
        if ii > 0 and ii % 100 == 0:
            sys.stdout.write('.')
            if ii % 1000 == 0:
                sys.stdout.write('%i/%i - kept %i\n' % (ii, sp.sum(is_candidate), len(keep_idx)))
            sys.stdout.flush()
        exon = ':'.join([data[i].chr, str(data[i].exons2[1, 0]), str(data[i].exons2[1, 1])]) 
        if not exon in gtex_exons[data[i].chr]:
            keep_idx.append(i)
            continue
        intron = ':'.join([data[i].chr, str(data[i].exons2[0, 1]), str(data[i].exons2[1, 0])])
        if intron in gtex_introns[data[i].chr]:
            continue
        intron = ':'.join([data[i].chr, str(data[i].exons2[1, 1]), str(data[i].exons2[2, 0])])
        if intron in gtex_introns[data[i].chr]:
            continue
        keep_idx.append(i)
    keep_idx = sp.array(keep_idx)
    print('%i of %i candidate events were not found in GTEx' % (keep_idx.shape[0], sp.sum(is_candidate)))
    is_candidate = sp.zeros(is_candidate.shape[0], dtype='bool')
    is_candidate[keep_idx] = 1
    pickle.dump(is_candidate, open(candidate_out_step2, 'wb'), -1)
else:
    print('%s already exists' % candidate_out_step2)
    is_candidate = pickle.load(open(candidate_out_step2, 'rb'), encoding='latin')
print('considering %i events' % sp.sum(is_candidate))

### recompute more stringent PSI for candidate events
if not os.path.exists(recomputed_psi):
    psis = []
    print('recomputing stricter PSI values')
    for ii,i in enumerate(sp.where(is_candidate)[0]):
        if ii > 0 and ii % 100 == 0:
            sys.stdout.write('.')
            if ii % 1000 == 0:
                sys.stdout.write('%i/%i\n' % (ii, sp.sum(is_candidate)))
            sys.stdout.flush()
        tmp = EV['event_counts'][:, :, i]
        a = tmp[:, [4, 5]].min(axis=1) 
        b = tmp[:, 6]
        nidx = sp.where(a + b < 10)[0]
        psis.append(a / (a + b))
        psis[-1][nidx] = sp.nan
    psis = sp.array(psis)
    pickle.dump(psis, open(recomputed_psi, 'wb'), -1)
else:
    print('load recomputed PSI values from %s' % recomputed_psi)
    psis = pickle.load(open(recomputed_psi, 'rb'), encoding='latin')

### filter by delta PSI
### collect other filter criteria such as
### - avg expression of the new exon
### - length of the new exon
### - total number of novel exon events in the same gene
psi_thresh = 0.1
keep_idx = []
cov_mean = []
cov_max = []
delta_psi = []
affected = []
non_affected = []
affected_files = []
non_affected_files = []
affected_donors = []
non_affected_donors = []
for ii,i in enumerate(sp.where(is_candidate)[0]):
    if ii > 0 and ii % 100 == 0:
        sys.stdout.write('.')
        if ii % 1000 == 0:
            sys.stdout.write('%i/%i - kept %i\n' % (ii, sp.sum(is_candidate), len(keep_idx)))
    sys.stdout.flush()
    #psi = EV['psi'][:, i][st_idx]
    psi = psis[ii, st_idx]
    nn_idx = sp.where(~sp.isnan(psi))[0]
    if nn_idx.shape[0] == 0:
        continue
    ### get number of samples that meet the psi threshold
    afct = (psi[nn_idx] > psi_thresh)
    if sp.sum(afct) == 0  :
        continue
    non_afct = (psi[nn_idx] <= psi_thresh)
    if sp.sum(non_afct) < 100:
        continue
    affected.append(sp.sum(afct))
    affected_files.append(','.join(files_ev[st_idx][nn_idx][afct][:10]))
    affected_donors.append(','.join(donors_ev[st_idx][nn_idx][afct][:10]))
    non_affected.append(sp.sum(non_afct))
    non_affected_files.append(','.join(files_ev[st_idx][nn_idx][non_afct][:10]))
    non_affected_donors.append(','.join(donors_ev[st_idx][nn_idx][non_afct][:10]))
    delta_psi.append(psi[nn_idx].max())
    keep_idx.append(i)
    cov_mean.append(sp.mean(EV['event_counts'][:, :, i][st_idx[nn_idx], 1]))
    cov_max.append(EV['event_counts'][:, :, i][st_idx[nn_idx], 1].max())
keep_idx = sp.array(keep_idx)
cov_mean = sp.array(cov_mean)
cov_max = sp.array(cov_max)
delta_psi = sp.array(delta_psi)
affected = sp.array(affected)
non_affected = sp.array(non_affected)
affected_files = sp.array(affected_files)
non_affected_files = sp.array(non_affected_files)
affected_donors = sp.array(affected_donors)
non_affected_donors = sp.array(non_affected_donors)

print('%i of %i candidate events did meet the delta psi threshold' % (keep_idx.shape[0], sp.sum(is_candidate)))
is_candidate = sp.zeros(is_candidate.shape[0], dtype='bool')
is_candidate[keep_idx] = 1
pickle.dump(is_candidate, open(candidate_out_step3, 'wb'), -1)

### report found candidates without checking for nearby variants
event_pos = sp.array([EV['event_pos'][i, :] for i in keep_idx])
gene_idx = sp.array([EV['gene_idx'][i] for i in keep_idx], dtype='int')
gene_strand = sp.array([EV['gene_strand'][i] for i in gene_idx])
gene_ids = sp.array([EV['gene_names'][i] for i in gene_idx])
gene_names = sp.array([un.get_ID(_, lookup=lookup) for _ in gene_ids])
is_coding = sp.array(['coding' if _ in coding_genes else 'non-coding' for _ in gene_ids])
event_chr = sp.array([EV['gene_chr'][i] for i in gene_idx])

### number of candidates per gene
gene_mult_dict = dict(sp.dstack(sp.unique(gene_ids, return_counts=True))[0, :])
gene_mult = sp.array([gene_mult_dict[_] for _ in gene_ids])

s_idx = sp.argsort(delta_psi)[::-1]
#           0         1          2--7       8         9           10         11           12         13        14       15        16            17         18              19
res = sp.c_[keep_idx, event_chr, event_pos, gene_ids, gene_names, delta_psi, gene_strand, is_coding, cov_mean, cov_max, affected, non_affected, gene_mult, affected_files, non_affected_files, affected_donors, non_affected_donors][s_idx, :]
sp.savetxt(os.path.join(BASEDIR_AS, 'alternative_splicing', 'exonization_candidates_C2.txt'), res, fmt='%s', delimiter='\t')


### continue overlapping to variants in a nearby window
keep_idx = []
### load genotype data
print('building variant index', file=sys.stderr)
trees = dict()
IN = h5py.File(VARIANTS, 'r')
pos_gt = IN['pos'][:]
strains_gt = IN['strains'][:]
chrms = sp.unique(pos_gt[:, 0])
for chrm in chrms:
    trees[chrm] = it.IntervalTree()
    idx = sp.where(pos_gt[:, 0] == chrm)[0]

    for i in range(0, idx.shape[0], 10000):
        begin = i
        end = min(i + 10000, idx.shape[0]) - 1
        trees[chrm][pos_gt[idx[begin], 1]:pos_gt[idx[end], 2] + 1] = idx[begin:end+1] 
print('loading sparse gt matrix', file=sys.stderr)
col = IN['gt_col'][:]
row = IN['gt_row'][:]
shape = IN['gt_shape'][:]
gt_in = spsp.coo_matrix((sp.array([True for _ in range(col.shape[0])], dtype='bool'), (row, col)), shape=shape).tocsc()
del col, row, shape
IN.close()

gt_keep = sp.where(sp.in1d(strains_gt, strains_ev[st_idx]))[0]
gt_in = gt_in[gt_keep, :]
strains_gt = strains_gt[gt_keep]

### do some strain matching
a, b = sp.where(strains_gt[:, sp.newaxis] == strains_ev[st_idx])
assert sp.all(strains_gt[a] == strains_ev[st_idx][b])

projects = sp.array([project_dict[_] for _ in strains_gt[a]])

keep_idx = []
delta_psi = []
alt_donors = []
ref_donors = []
alt_projects = []
var_pos = []
alt_files = []
ref_files = []

for ii,i in enumerate(sp.where(is_candidate)[0]):
    if ii > 0 and ii % 100 == 0:
        sys.stdout.write('.')
        if ii % 1000 == 0:
            sys.stdout.write('%i/%i - kept %i\n' % (ii, sp.sum(is_candidate), len(keep_idx)))
        sys.stdout.flush()

    ### get variant and alt strains
    chrm = data[i].chr
    if chrm == 'X':
        chrm = '23'
    if chrm == 'Y':
        chrm = '24'
    pos11 = data[i].exons2[1, 0] - 25
    pos12 = data[i].exons2[1, 0] + 25
    pos21 = data[i].exons2[1, 1] - 25
    pos22 = data[i].exons2[1, 1] + 25
    idx = []
    for itv in trees[int(chrm)].search(pos11, pos12):
        _idx = sp.where((pos11 < pos_gt[itv.data, 1]) & (pos12 > pos_gt[itv.data, 2]))[0]
        if _idx.shape[0] > 0:
            idx.append(itv.data[_idx])
    for itv in trees[int(chrm)].search(pos21, pos22):
        _idx = sp.where((pos21 < pos_gt[itv.data, 1]) & (pos22 > pos_gt[itv.data, 2]))[0]
        if _idx.shape[0] > 0:
            idx.append(itv.data[_idx])
    idx = sp.unique([x for _ in idx for x in _])
    if len(idx) < 1:
        continue
    for jdx in idx:
        gt = gt_in[:, jdx].toarray()[a, :]
        if sp.sum(gt) == 0:
            continue
        ### check whether we can get a readout from PSI values at the variant positions
        alt_idx = sp.where(gt)[0]
        ref_idx = sp.where(~gt)[0]
        if sp.all(sp.isnan(EV['psi'][:, i][st_idx][b][alt_idx])):
            continue
        if EV['event_counts'][:, :, i][st_idx[b][alt_idx], :][:, [4, 5]].min() < read_thresh:
            continue
        keep_idx.append(i)
        delta_psi.append(sp.absolute(spst.nanmean(EV['psi'][:, i][st_idx][b][alt_idx]) - spst.nanmean(EV['psi'][:, i][st_idx][b][ref_idx])))
        alt_donors.append(','.join(strains_gt[a][alt_idx]))
        pro = []
        ref = []
        ref_donor = []
        alt = []
        for _ in alt_idx:
            pro.append(projects[_])
            iidx = sp.where(projects == projects[_])[0]
            iidx = iidx[~sp.in1d(iidx, alt_idx)]
            alt.append(files_ev[st_idx][b][_])
            try:
                ref.append(files_ev[st_idx][b][iidx][0])
                ref_donor.append(strains_gt[a][iidx][0])
            except:
                ref.append('NA')
                ref_donor.append('NA')
        alt_projects.append(','.join(pro))
        var_pos.append(':'.join(pos_gt[jdx, :].astype('str')))
        alt_files.append(','.join(alt))
        ref_files.append(','.join(ref))
        ref_donors.append(','.join(ref_donor))

keep_idx = sp.array(keep_idx)
print('%i of %i candidate events are in context of a somatic alteration' % (keep_idx.shape[0], sp.sum(is_candidate)))
is_candidate = sp.zeros(is_candidate.shape[0], dtype='bool')
is_candidate[keep_idx] = 1
pickle.dump(is_candidate, open(candidate_out_step4, 'wb'), -1)

delta_psi = sp.array(delta_psi)
alt_donors = sp.array(alt_donors)
ref_donors = sp.array(ref_donors)
alt_projects = sp.array(alt_projects)
var_pos = sp.array(var_pos)
alt_files = sp.array(alt_files)
ref_files = sp.array(ref_files)

k_idx = sp.where(sp.in1d(keep_idx, conf_idx))[0]
keep_idx = keep_idx[k_idx]
delta_psi = delta_psi[k_idx]
alt_donors = alt_donors[k_idx]
ref_donors = ref_donors[k_idx]
alt_projects = alt_projects[k_idx]
var_pos = var_pos[k_idx]
alt_files = alt_files[k_idx]
ref_files = ref_files[k_idx]

event_pos = sp.array([EV['event_pos'][i, :] for i in keep_idx])
gene_idx = sp.array([EV['gene_idx'][i] for i in keep_idx], dtype='int')
gene_strand = sp.array([EV['gene_strand'][i] for i in gene_idx])
gene_ids = sp.array([EV['gene_names'][i] for i in gene_idx])
gene_names = sp.array([un.get_ID(_, lookup=lookup) for _ in gene_ids])
is_coding = sp.array(['coding' if _ in coding_genes else 'non-coding' for _ in gene_ids])
event_chr = sp.array([EV['gene_chr'][i] for i in gene_idx])

s_idx = sp.argsort(delta_psi)[::-1]
#           0        1         2          3--8       9         10          11         12          13          14            15         16         17           18
res = sp.c_[var_pos, keep_idx, event_chr, event_pos, gene_ids, gene_names, delta_psi, alt_donors, ref_donors, alt_projects, alt_files, ref_files, gene_strand, is_coding][s_idx, :]
sp.savetxt(os.path.join(BASEDIR_AS, 'alternative_splicing', 'exonization_candidates_C2.SVovlp.txt'), res, fmt='%s', delimiter='\t')
EV.close()


