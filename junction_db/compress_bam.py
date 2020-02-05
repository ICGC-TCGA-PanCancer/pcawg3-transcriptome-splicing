import sys
import h5py
import re
import pysam
import scipy as sp
import scipy.sparse

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'MANDATORY')
    required.add_option('-b', '--bams', dest='bams', metavar='FILE1,FILE2,...', help='alignment files in BAM format (comma separated list)', default='-')
    output = OptionGroup(parser, 'OPTIONS')
    output.add_option('-v', '--verbose', dest='verbose', metavar='y|n', help='verbosity', default='n')
    output.add_option('-p', '--parallel', dest='parallel', metavar='<INT>', type='int', help='use multiple processors [1]', default=1)
    parser.add_option_group(required)
    parser.add_option_group(output)

    (options, args) = parser.parse_args()

    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    options.parser = parser
    return options


def parse_header(header_string):

    hd = dict()

    for line in header_string.strip('\n').split('\n'):
        sl = line.strip().split('\t')
        ### ignore comment lines
        if sl[0] == '@CO':
            continue
        td = dict([x.split(':', 1) for x in sl[1:] if ':' in x])    
        try:
            hd[sl[0].strip('@')].append(td)
        except KeyError:
            hd[sl[0].strip('@')] = [td]
    return hd


def sort_rows(array, index = None):
    """Sort array by rows"""

    ### empty array
    if array.shape[0] == 0:
        if index == True:
            return (array, [])
        else:
            return (array)

    ### only one row
    if len(array.shape) == 1:
        if index == True:
            return (array, [0])
        else:
            return (array)

    ### more than one row
    s_idx = sp.lexsort([array[:, -i] for i in range(1, array.shape[1] + 1)])

    if index == True:
        return (array[s_idx, :], s_idx)
    else:
        return array[s_idx, :]


def filter_read(read, filt, spliced, mapped, strand, primary_only, no_mm=False):

    if read.is_unmapped:
        return True
    if primary_only and read.is_secondary:
        return True

    is_spliced = ('N' in read.cigarstring)
    if is_spliced:
        if not spliced:
            return True
    elif not mapped:
        return True

    tags = dict(read.tags)

    if filt is not None:
        ### handle mismatches
        if not no_mm:
            try:
                if filt['mismatch'] < tags['NM']:
                    return True
            except KeyError:
                print('Expect the NM tag to be present for filtering', file=sys.stderr)
                sys.exit(1)

        if is_spliced:
            ### handle min segment length
            ### remove all elements from CIGAR sting that do not 
            ### contribute to the segments (hard- and softclips and insertions)
            cig = re.sub(r'[0-9]*[HSI]', '', read.cigarstring)
            ### split the string at the introns and sum the remaining segment elements, compare to filter
            if min([sum([int(y) for y in re.split('[MD]', x)[:-1]]) for x in re.split('[0-9]*N', cig)]) <= filt['exon_len']:
                return True
    
    ### check strand information
    if strand is not None:
        try:
            if tags['XS'] != strand:
                return True
        except KeyError:
            pass

    return False


def summarize_chr(fname, chr_name, verbose=False, filt=None, strand=None, mapped=True, spliced=True, unstranded=True, primary_only=True):

    infile = pysam.Samfile(fname, 'rb')
    introns_p = dict()
    introns_m = dict()

    if verbose:
        print('Summarizing contig %s of file %s' % (chr_name, fname), file=sys.stdout)

    chr_len = [int(x['LN']) for x in parse_header(infile.text)['SQ'] if x['SN'] == chr_name]
    if len(chr_len) == 0:
        print('No information found for contig %s' % (chr_name), file=sys.stdout)
        return (chr_name, scipy.sparse.coo_matrix(sp.zeros((0, 1)), dtype='uint32'), sp.zeros((0, 3), dtype='uint32'), sp.zeros((0, 3), dtype='uint32'))
    chr_len = chr_len[0]

    ### read matrix has three rows: 0 - no strand info, 1 - plus strand info, 2 - minus strand info
    read_matrix = sp.zeros((3, chr_len), dtype='uint32') 

    if infile.gettid(chr_name) > -1:
        ### pysam query is zero based in position (results are as well), all intervals are pythonic half open
        for read in infile.fetch(chr_name, until_eof=True):

            ### check if we skip this reads
            if filter_read(read, filt, spliced, mapped, strand, primary_only):
                continue
            
            tags = dict(read.tags)
            curr_read_stranded = ('XS' in tags)
            is_minus = False
            if curr_read_stranded:
                is_minus = (tags['XS'] == '-')
                    
            ### get introns and covergae
            p = read.pos 
            for o in read.cigar:
                if o[0] == 3:
                    if is_minus:
                        try:
                            introns_m[(p, p + o[1])] += 1
                        except KeyError:
                            introns_m[(p, p + o[1])] = 1
                    else:
                        try:
                            introns_p[(p, p + o[1])] += 1
                        except KeyError:
                            introns_p[(p, p + o[1])] = 1
                if o[0] in [0, 2]:
                    read_matrix[0 + int(curr_read_stranded) + int(is_minus), p:(p + o[1])] += 1
                if o[0] in [0, 2, 3]:
                    p += o[1]

    ### convert introns into scipy array
    if len(introns_p) >= 1:
        introns_p = sp.array([[k[0], k[1], v] for k, v in introns_p.items()], dtype='uint32')
        introns_p = sort_rows(introns_p)
    else:
        introns_p = sp.zeros(shape=(0, 3), dtype='uint32')
    if len(introns_m) >= 1:
        introns_m = sp.array([[k[0], k[1], v] for k, v in introns_m.items()], dtype='uint32')
        introns_m = sort_rows(introns_m)
    else:
        introns_m = sp.zeros(shape=(0, 3), dtype='uint32')

    ### make read matrix sparse
    if unstranded:
        read_matrix = scipy.sparse.coo_matrix(sp.sum(read_matrix, axis=0)[sp.newaxis, :], dtype='uint32')
    else:
        read_matrix = scipy.sparse.coo_matrix(read_matrix, dtype='uint32')

    return (chr_name, read_matrix, introns_m, introns_p)



def main():

    ### get command line options
    options = parse_options(sys.argv)
    if options.bams == '-':
        options.parser.print_help()
        sys.exit()
        

    ### build read filter dict
    read_filter = None

    for bfn in options.bams.split(','):

        if options.verbose == 'y':
            print('processing %s' % bfn)

        infile = pysam.Samfile(bfn, 'rb')
        chrms = [x['SN'] for x in parse_header(infile.text)['SQ']]
        infile.close()

        OUT = h5py.File(re.sub(r'.bam$', '', bfn) + '.hdf5', 'w')
        if options.parallel > 1:
            import multiprocessing as mp
            pool = mp.Pool(processes=options.parallel)
            result = [pool.apply_async(summarize_chr, args=(bfn, str(chrm),), kwds={'filt':read_filter}) for chrm in sorted(chrms)]
            while result:
                tmp = result.pop(0).get()
                if options.verbose == 'y':
                    print(tmp[0])
                OUT.create_dataset(name=(tmp[0] + '_reads_row'), data=tmp[1].row.astype('uint8'), compression='gzip')
                OUT.create_dataset(name=(tmp[0] + '_reads_col'), data=tmp[1].col, compression='gzip')
                OUT.create_dataset(name=(tmp[0] + '_reads_dat'), data=tmp[1].data, compression='gzip')
                OUT.create_dataset(name=(tmp[0] + '_reads_shp'), data=tmp[1].shape)
                OUT.create_dataset(name=(tmp[0] + '_introns_m'), data=tmp[2], compression='gzip')
                OUT.create_dataset(name=(tmp[0] + '_introns_p'), data=tmp[3], compression='gzip')
                del tmp
        else:
            for chrm in sorted(chrms):
                if options.verbose == 'y':
                    print(chrm)
                tmp = summarize_chr(bfn, str(chrm), filt=read_filter)
                OUT.create_dataset(name=(chrm + '_reads_row'), data=tmp[1].row.astype('uint8'), compression='gzip')
                OUT.create_dataset(name=(chrm + '_reads_col'), data=tmp[1].col, compression='gzip')
                OUT.create_dataset(name=(chrm + '_reads_dat'), data=tmp[1].data, compression='gzip')
                OUT.create_dataset(name=(chrm + '_reads_shp'), data=tmp[1].shape)
                OUT.create_dataset(name=(chrm + '_introns_m'), data=tmp[2], compression='gzip')
                OUT.create_dataset(name=(chrm + '_introns_p'), data=tmp[3], compression='gzip')
        OUT.close()

if __name__ == "__main__":
    main()
