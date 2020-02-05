import sys
import os
import scipy as sp

sys.path.append('..')
from paths import BASEDIR_AS
basedir = os.path.join(BASEDIR_AS, 'alternative_splicing')

full_data = sp.loadtxt(os.path.join(basedir, 'exonization_candidates_C2.txt'), dtype='str', delimiter='\t')
snv_data = sp.loadtxt(os.path.join(basedir, 'exonization_candidates_C2.SVovlp.txt'), dtype='str', delimiter='\t')

fin_data = []
### keep relevant subset of data
### event id
fin_data.append(full_data[:, 0])
### event pos
fin_data.append(sp.array([x[1] + '-' + ':'.join(x[2:8])  for x in full_data]))
### strand
fin_data.append(full_data[:, 11])
### ensemble id
fin_data.append(full_data[:, 8])
### gene name
fin_data.append(full_data[:, 9])
### max dPSI
fin_data.append(full_data[:, 10])
### coding status
fin_data.append(full_data[:, 12])
### overlapping SNVs
snvs = []
for i in range(full_data.shape[0]):
    tmp = []
    for j in sp.where(snv_data[:, 1] == full_data[i, 0])[0]:
        tmp.append(snv_data[j, 0])
    if len(tmp) > 0:
        snvs.append(','.join(tmp))
    else:
        snvs.append('NA')
fin_data.append(sp.array(snvs))

### gen header
header = sp.array(['event_id', 'event_pos', 'strand', 'ensemble_id', 'gene_name', 'max_dpsi', 'coding_status', 'overlap_snv'])
fin_data = sp.r_[header[sp.newaxis, :], sp.array(fin_data).T]
sp.savetxt(os.path.join(basedir, 'supplemental_table_exonization_candidates.tsv'), fin_data, fmt='%s', delimiter='\t') 
