#! /usr/bin/env python 
import sys
import os
import re
import h5py
import scipy as sp
import pickle

from modules.hdf5 import appendToHDF5

#event_types = ['exon_skip', 'intron_retention', 'alt_3prime', 'alt_5prime', 'mutex_exons', 'mult_exon_skip']

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'MANDATORY')
    required.add_option('-i', '--indirs', dest='indirs', metavar='DIR1,DIR2,...', help='comma separated list of spladder directories to be merged', default='-')
    required.add_option('-o', '--outdir', dest='outdir', metavar='DIR', help='spladder directory containing the integrated results of the input directories', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-c', '--confidence', dest='confidence', metavar='INT', type='int', help='confidence level (0 lowest to 3 highest) [3]', default=3)
    optional.add_option('-t', '--event_types', dest='event_types', metavar='TYPE1,TYPE2,...', help='list of parts considered for merging [genes,exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons]', default='genes,exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons')
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default=False)

    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()
    options.event_types = options.event_types.strip(',').split(',')

    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    options.parser = parser
    return options


def main():
    
    """Workflow to merge two SplAdder directories into a new one"""
    
    ### parse command line parameters
    options = parse_options(sys.argv)
    
    ### create plot directory if it does not exist yet
    if not os.path.exists(options.outdir):
        os.mkdir(options.outdir)
    else:
        print('WARNING: %s already exists - any content will be overwritten!' % options.outdir, file=sys.stderr)

    ### loop over input directories
    for i, indir in enumerate(options.indirs.strip(',').split(',')):
        if options.verbose:
            print('\nhandling directory %s' % indir, file=sys.stderr)

        if 'genes' in options.event_types:
            ### handle merged gene files
            infile = os.path.join(indir, 'spladder', 'genes_graph_conf%i.merge_graphs.count.hdf5' % (options.confidence))
            if options.verbose:
                print('processing %s' % infile, file=sys.stderr)
            IN = h5py.File(infile, 'r')
            outdir = os.path.join(options.outdir, 'spladder')
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            OUT = h5py.File(os.path.join(outdir, 'genes_graph_conf%i.merge_graphs.count.hdf5' % (options.confidence)), 'a')
            if i == 0:
                for key in list(IN.keys()):
                    if options.verbose:
                        print('\tkey: %s' % key, file=sys.stderr)
                    if key in ['edge_idx', 'gene_ids_edges', 'gene_ids_segs', 'gene_names', 'seg_len']:
                        OUT.create_dataset(name=key, data=IN[key][:], compression='gzip', chunks=True) 
                    elif key in ['strains']:
                        OUT.create_dataset(name=key, data=IN[key][:], compression='gzip', chunks=True, maxshape=(None,)) 
                    else:
                        OUT.create_dataset(name=key, data=IN[key][:], compression='gzip', chunks=True, maxshape=(IN[key].shape[0], None)) 
            else:
                for key in list(IN.keys()):
                    if options.verbose:
                        print('\tkey: %s' % key, file=sys.stderr)
                    if key in ['edge_idx', 'gene_ids_edges', 'gene_ids_segs', 'gene_names', 'seg_len']:
                        assert IN[key].shape[0] == OUT[key].shape[0]
                    elif key in ['strains']:
                        appendToHDF5(OUT, IN[key][:], key)
                    else:
                        appendToHDF5(OUT, IN[key][:], key, faxis=1, daxis=1)
            OUT.close()
            IN.close()
                    
        ### handle event files
        for event_type in options.event_types:
            if event_type == 'genes':
                continue
            infile = os.path.join(indir, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, options.confidence))
            if options.verbose:
                print('processing %s' % infile, file=sys.stderr)
            if not os.path.exists(infile):
                continue
            IN = h5py.File(infile, 'r')
            OUT = h5py.File(os.path.join(options.outdir, os.path.basename(infile)), 'a')
            for key in list(IN.keys()):
                if options.verbose:
                    print('\tkey: %s' % key, file=sys.stderr)
                if not key in OUT:
                    if key == 'event_features':
                        GRP = OUT.create_group(name=key)
                        for kkey in list(IN[key].keys()):
                            GRP.create_dataset(name=kkey, data=IN[key][kkey][:])
                    elif key == 'psi':
                        OUT.create_dataset(name=key, data=IN[key][:], compression='gzip', chunks=True, maxshape=(None, IN[key].shape[1]))
                    elif key == 'verified':
                        OUT.create_dataset(name=key, data=IN[key][:], compression='gzip', chunks=True, maxshape=(IN[key].shape[0], None, IN[key].shape[2]))
                    elif key == 'strains':
                        OUT.create_dataset(name=key, data=IN[key][:], compression='gzip', chunks=True, maxshape=(None,))
                    elif key == 'event_counts':
                        if len(IN[key].shape) == 1:
                            continue
                        else:
                            chunksize = IN[key].chunks[0]
                            for c, chunk in enumerate(range(0, IN[key].shape[0], chunksize)):
                                if options.verbose:
                                    sys.stdout.write('%i/%i chunks done\n' % (c, IN[key].shape[0] / chunksize + 1))
                                    sys.stdout.flush()
                                curr_range = list(range(chunk, min(chunk + chunksize, IN[key].shape[0])))
                                if c == 0:
                                    OUT.create_dataset(name=key, data=IN[key][curr_range, :, :], compression='gzip', chunks=True, maxshape=(None, IN[key].shape[1], IN[key].shape[2]))
                                else:
                                    appendToHDF5(OUT, IN[key][curr_range, :, :], key, faxis=0, daxis=0) 
                    else:
                        OUT.create_dataset(name=key, data=IN[key][:], compression='gzip', chunks=True)
                else:
                    if key == 'event_features':
                        continue
                    elif key in ['psi', 'strains']:
                        appendToHDF5(OUT, IN[key][:],key)
                    elif key == 'verified':
                        appendToHDF5(OUT, IN[key][:], key, faxis=1, daxis=1)
                    elif key == 'event_counts':
                        if len(IN[key].shape) == 1:
                            continue
                        else:
                            chunksize = IN[key].chunks[0]
                            for c, chunk in enumerate(range(0, IN[key].shape[0], chunksize)):
                                if options.verbose:
                                    sys.stdout.write('%i/%i chunks done\n' % (c, IN[key].shape[0] / chunksize + 1))
                                    sys.stdout.flush()
                                curr_range = list(range(chunk, min(chunk + chunksize, IN[key].shape[0])))
                                appendToHDF5(OUT, IN[key][curr_range, :, :], key, faxis=0, daxis=0) 
                    elif i > 0 and key in ['conf_idx', 'confirmed', 'num_verified']:
                        continue
                    else:
                        assert sp.all(IN[key].shape == OUT[key].shape)

                ### handle confirmed keys
                if i > 0 and 'verified' in list(IN.keys()):
                    ### num_verified
                    vshape = OUT['verified'].shape
                    num_verified = sp.empty((vshape[0], vshape[2]), dtype='int')
                    for j in range(0, vshape[0], 1000):
                        idx = list(range(j, min(vshape[0], j + 1000)))
                        num_verified[idx, :] = OUT['verified'][idx, :, :].sum(axis=1)
                    del OUT['num_verified']
                    OUT.create_dataset(name='num_verified', data=num_verified, compression='gzip')

                    ### confirmed
                    del OUT['confirmed']
                    confirmed = num_verified.min(axis=1)
                    OUT.create_dataset(name='confirmed', data=confirmed, compression='gzip')

                    ### confidx
                    if 'conf_idx' in list(IN.keys()):
                        conf_idx = sp.where(confirmed > 0)[0]
                        del OUT['conf_idx']
                        OUT.create_dataset(name='conf_idx', data=conf_idx, compression='gzip')
                        del conf_idx

                    ### clean up
                    del num_verified
                    del confirmed
            OUT.close()
            IN.close()


if __name__ == "__main__":
    main()

