import sys
import os
import pickle
import h5py
import scipy as sp
import re
import glob

if len(sys.argv) < 5:
    sys.stderr.write('Usage: %s <coordinate_db> <projected_dir> <threshold> <outfile> \n' % sys.argv[0])
    sys.exit(1)

junc_pickle = sys.argv[1]
projected_dir = sys.argv[2]
thresh = sys.argv[3]
outfname = sys.argv[4]

### load coordinate DB
junctions = pickle.load(open(junc_pickle, 'rb'))

### collect files containing projected counts
flist = glob.glob(os.path.join(projected_dir, '*.t%s.projected.hdf5' % thresh))

### gather single result files into joint DB
OUT = h5py.File(outfname, 'w')
strains = []
for i,fname in enumerate(flist):
    sys.stderr.write('processing %i/%i\n' % (i + 1, len(flist)))
    strains.append(os.path.basename(fname).split('.')[0])
    with h5py.File(fname, 'r') as IN:
        for key in IN:
            if IN[key].shape[0] == 0:
                continue
            if not key in OUT:
                OUT.create_dataset(name=key, data=IN[key][:], compression='gzip', maxshape=(IN[key].shape[0], None)) 
                ### add positions
                if key.endswith('_m'):
                    dbkey = (key.split('_')[1], 2)
                elif key.endswith('_p'):
                    dbkey = (key.split('_')[1], 1)
                pos = sp.array([[_[0], _[1], junctions[dbkey][_]]for _ in junctions[dbkey]])
                sidx = sp.argsort(pos[:, 2])
                pos = pos[sidx, :]
                OUT.create_dataset(name=re.sub(r'^junctions_', 'pos_', key), data=sp.array(pos)[:, :2], compression='gzip')
            else:
                tmp = OUT[key].shape
                OUT[key].resize((tmp[0], tmp[1] + 1))
                OUT[key][:, tmp[1]] = IN[key][:, 0]
        ### generate entries for all contigs not present in the file
        for dbkey in junctions:
            if len(junctions[dbkey]) == 0:
                continue
            if dbkey[1] == 2:
                key = 'junctions_%s_m' % dbkey[0]
            else:
                key = 'junctions_%s_p' % dbkey[0]
            if not key in IN:
                if not key in OUT:
                    OUT.create_dataset(name=key, data=sp.zeros((len(junctions[dbkey]), 1), dtype='int'), compression='gzip', maxshape=(len(junctions[dbkey]), None))
                    pos = sp.array([[_[0], _[1], junctions[dbkey][_]]for _ in junctions[dbkey]])
                    sidx = sp.argsort(pos[:, 2])
                    pos = pos[sidx, :]
                    OUT.create_dataset(name=re.sub(r'^junctions_', 'pos_', key), data=sp.array(pos)[:, :2], compression='gzip')
                else:
                    tmp = OUT[key].shape
                    OUT[key].resize((tmp[0], tmp[1] + 1))
                    OUT[key][:, tmp[1]] = sp.zeros((len(junctions[dbkey]),), dtype='int')
                    
OUT.create_dataset(name='strains', data=sp.array(strains, dtype='str').view(sp.chararray).encode('utf-8'), compression='gzip')
OUT.close()

