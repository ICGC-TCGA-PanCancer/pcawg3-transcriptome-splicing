import scipy as sp
import pickle
import os
import pdb


def get_tags_gtf(tagline):
    """Extract tags from given tagline"""

    tags = dict()
    for t in tagline.strip(';').split(';'):
        tt = t.strip(' ').split(' ')
        tags[tt[0]] = tt[1].strip('"')
    return tags


def parse_gtf(infile):

    ### init dicts
    lookup = dict()
    lookup['ensembl2genename'] = dict()
    lookup['ensembl2longensembl'] = dict()
    lookup['ensembl2genetype'] = dict()
    lookup['ensembl2havana'] = dict()
    lookup['ensembl2transcript'] = dict()
    lookup['genename2ensembl'] = dict()
    lookup['transcript2gene'] = dict()

    print('Parsing annotation data from %s (may take a while) ...' % infile)

    ### parse input gtf file and gather information
    for line in open(infile, 'r'):
        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        
        if not sl[2] in ['gene', 'transcript']:
            continue
        
        tags = get_tags_gtf(sl[8])

        gene_id = tags['gene_id'].split('.')[0]
        if sl[2] == 'gene':
            lookup['ensembl2genename'][gene_id] = tags['gene_name']
            lookup['ensembl2longensembl'][gene_id] = tags['gene_id']
            lookup['ensembl2genetype'][gene_id] = tags['gene_type']
            try:
                lookup['ensembl2havana'][gene_id] = tags['havana_gene']
            except KeyError:
                lookup['ensembl2havana'][gene_id] = 'NA'
            lookup['genename2ensembl'][tags['gene_name']] = gene_id
        if sl[2] == 'transcript':
            try:
                lookup['ensembl2transcript'][gene_id].append(tags['transcript_id'])
            except KeyError:
                lookup['ensembl2transcript'][gene_id] = [tags['transcript_id']]
            lookup['transcript2gene'][tags['transcript_id']] = gene_id
    print('...done')

    return lookup


def get_lookup_complete():

    basedir = '/cluster/work/grlab/projects/ICGC/annotation'

    ### get lookup tables from annotation
    if os.path.exists('%s/gencode.v19.annotation.hs37d5_chr.gtf.lookup.pickle' % basedir):
        lookup = pickle.load(open('%s/gencode.v19.annotation.hs37d5_chr.gtf.lookup.pickle' % basedir, 'rb'), encoding='latin')
    else:
        lookup = parse_gtf(os.path.join(basedir, 'gencode.v19.annotation.hs37d5_chr.gtf'))
        pickle.dump(lookup, open('%s/gencode.v19.annotation.hs37d5_chr.gtf.lookup.pickle' % basedir, 'wb'), -1)

    ### get lookup table from hgnc mapping
    data = sp.loadtxt('%s/gencode.v19.annotation.ensembl2hgnc.tsv' % basedir, delimiter='\t', dtype='str')
    lookup['ensembl2hgnc'] = dict()
    for i in range(data.shape[0]):
        #try:
        #    lookup['ensembl2hgnc'][data[i, 0]].append(data[i, 2] if len(data[i, 2]) > 0 else 'NA')
        #except KeyError:
        lookup['ensembl2hgnc'][data[i, 0]] = [data[i, 2] if len(data[i, 2]) > 0 else 'NA']

    return lookup


def get_ID(en_id, lookup=None, which='name', aliases=None):
    
    if lookup is None:
        lookup = get_lookup_complete()

    if which != 'ensembl':
        en_id_ = en_id.split('.')[0]
    else:
        en_id_ = en_id

    ret = en_id
    alias_idx = -1

    while ret == en_id:
        #print 'looking up %s' % en_id
        if (alias_idx >= 0):
            if alias_idx >= len(alias):
                break
            en_id_ = alias[alias_idx]
            alias_idx += 1

        try:
            if which == 'name':
                ret = lookup['ensembl2genename'][en_id_]
            elif which == 'havana':
                ret = lookup['ensembl2havana'][en_id_]
            elif which == 'hgnc':
                ret = lookup['ensembl2hgnc'][en_id_]
            elif which == 'type':
                ret = lookup['ensembl2genetype'][en_id_]
            elif which == 'transcripts':
                ret = lookup['ensembl2transcript'][en_id_]
            elif which == 'ensembl':
                ret = lookup['genename2ensembl'][en_id_]
            elif which == 'uniprot':
                ret = get_uniprotID(en_id)
        except:
            ### check whether we have any aliases for this gene alias
            if aliases is None:
                break
            if alias_idx < 0:
                idx = sp.where([en_id in x.split('|') for x in aliases])[0]
                if idx.shape[0] == 0:
                    break
                alias = aliases[idx[0]].split('|')
                alias_idx = 0

    return ret


if __name__ == "__main__":
    lookup = get_lookup_complete()
    print(get_ID('ENSG00000063244', lookup=lookup))
    print(get_ID('ENSG00000063244.21', lookup=lookup))
    print(get_ID('ENSG00000063244', lookup=lookup, which='type'))
    print(get_ID('ENSG00000063244', lookup=lookup, which='havana'))
    print(get_ID('ENSG00000063244', lookup=lookup, which='transcripts'))
    print(get_ID('ALK', which='ensembl', lookup=lookup))

