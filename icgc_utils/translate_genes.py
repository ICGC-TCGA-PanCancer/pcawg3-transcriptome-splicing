import sys
import os
import utils.names as un

if len(sys.argv) < 2:
    print >> sys.stdout, 'Usage: %s <genelist.txt>'
    sys.exit(1)

l = un.get_lookup_complete()
for line in open(sys.argv[1]):
    sl = line.strip().split('\t')
    print un.get_ID(sl[0], lookup=l)
