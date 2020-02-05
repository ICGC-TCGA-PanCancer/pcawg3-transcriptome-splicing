import sys
import scipy as sp
import gzip
import re

if len(sys.argv) < 2:
    print('Usage: %s <icgc.txt.gz - file>' % sys.argv[0], file=sys.stderr)
    sys.exit(1)
infname = sys.argv[1]

from paths import BASEDIR

### get metadata 
metatable = os.path.join(BASEDIR, 'orig_data/metadata/per_aliquot_v2/rnaseq_metadata.tsv')
metadata = []
for line in open(metatable, 'r'):
    sl = line.strip().split('\t')
    metadata.append(sl)
metadata = sp.array(metadata, dtype='str')
header = metadata[0, :]
metadata = metadata[1:, :]
aidx = sp.where(header == 'analysis_id')[0]
alidx = sp.where(header == 'aliquot_id')[0]

outfname = re.sub(r'txt.gz$', '', infname) + 'aliquotID.txt.gz'
out = gzip.open(outfname, 'w')

for i, line in enumerate(gzip.open(infname, 'r')):
    if i > 1 and i % 100 == 0:
        sys.stdout.write('.')
        if i % 1000 == 0:
            sys.stdout.write('%i\n' % i)
        sys.stdout.flush()
    if i == 0:
        header = line.strip().split('\t')
        strains = sp.array([x.split('.')[0] for x in header[6:]], dtype='str')
        offset = len(header) - len(strains)
        a,b = sp.where(strains[:, sp.newaxis] == metadata[:, aidx[0]])
        assert sp.all(a == sp.sort(a))
        assert sp.all(strains[a] == metadata[b, aidx[0]])
        header = sp.r_[header[:offset], metadata[b, alidx[0]]]
        k_idx = sp.r_[sp.arange(offset), a + offset]
        print('\t'.join(header), file=out)
    else:
        sl = sp.array(line.strip().split('\t'))
        print('\t'.join(sl[k_idx]), file=out)
out.close()
