import sys
import os
import scipy as sp
import scipy.sparse as spsp
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pickle
import re
import numpy.random as npr
npr.seed(23)
import intervaltree as it
sys.path.append('/cluster/project/grlab/home/akahles/git/software/spladder/python')
import modules.viz.coverage as cv
import modules.viz.genelets as gl
sys.path.append('/cluster/project/grlab/home/akahles/git/tools/python/viz/')
import axes as axs

sys.path.append('../icgc_anno')
import translate_metadata as tm

sys.path.append('../icgc_utils')
import names as un

if len(sys.argv) < 2:
   print('Usage: %s <result_file_exonization>' % sys.argv[0], file=sys.stderr)
   sys.exit(1)
fname = sys.argv[1]

TOP = 100
CONF = 2
PLOTDIR = '/cluster/project/grlab/projects/ICGC/alternative_splicing_new/plots/exonization_detail_noSV_C%s' % CONF
if not os.path.exists(PLOTDIR):
    os.makedirs(PLOTDIR)

results = sp.loadtxt(fname, dtype='str', delimiter='\t')
if len(results.shape) < 2:
    results = results[sp.newaxis, :]

lookup = un.get_lookup_complete()
logtag = ''
if ('--log' in sys.argv):
    logtag = '.log'

gs = gridspec.GridSpec(2, 2, width_ratios=[5, 1], height_ratios=[2, 1])

[strain_dict, tumor_dict, primary_dict, donor_dict] = tm.translate([('analysis_id', 'icgc_donor_id'), ('analysis_id', 'is_tumour'), ('analysis_id', 'specimen_type'), ('icgc_donor_id', 'analysis_id')])

for i in range(min(TOP, results.shape[0])):

    print('Processing %i: %s' % (i, results[i, 0]), file=sys.stderr)

    alt_donors = sp.array(results[i, 20].split(','))
    ref_donors = sp.array(results[i, 21].split(','))
    #projects = results[i, 14].split(',')
    files_alt = results[i, 18].split(',')
    files_ref = results[i, 19].split(',')

    #var_pos = [int(x.split(':')[1]) for x in results[i, 0].split(',')]

    ### generate intron filter
    event_pos = sp.array([int(results[i, 2]) - 100, int(results[i, 7]) + 100])
    exons2 = sp.reshape(results[i, 2:8], (3, 2)).astype('int')
    exons1 = exons2[[0, 2], :]
    if (exons2[1, 1] - exons2[1, 0]) % 3 == 0:
        in_frame = 'yes'
    else:
        in_frame = 'no'
    intron_filter = [(int(results[i, 3]), int(results[i, 4])),\
                     (int(results[i, 5]), int(results[i, 6])),\
                     (int(results[i, 3]), int(results[i, 6]))]
    event_chr = str(results[i, 1])

    rows = alt_donors.shape[0] + 1
    height_ratios = [3 for _ in range(rows)]
    height_ratios[-1] = 1
    gs = gridspec.GridSpec(rows, 1, height_ratios=height_ratios)
    fig = plt.figure(figsize = (16, 3 * rows), dpi=200)
    #fig.suptitle('%s in %s (Variant: %s:%s)' % (results[i, 6], results[i, 7],results[i, 0],  results[i, 2]), fontsize=12, y=1.03)
        
    for d, donor in enumerate(alt_donors):

        ### plot coverage profile
        ax = fig.add_subplot(gs[d, 0])
        cax = []
        labels = []
        #ax.set_title('%s - Gene: %s (%s) - SNV: %s - in frame: %s' % (results[i, 14], un.get_ID(results[i, 9], lookup=lookup), results[i, 17], results[i, 0], in_frame), fontsize=10)
        ax.set_title('Gene: %s -in frame: %s' % (un.get_ID(results[i, 8], lookup=lookup), in_frame), fontsize=10)
        cax.append(cv.cov_from_bam(event_chr, event_pos.min(), event_pos.max(), [files_alt[d]], color_cov='r', ax=ax, intron_cnt=True, color_intron_edge='r', log=('--log' in sys.argv), intron_filter=intron_filter, return_legend_handle=True, label='Alt: %s' % donor))
        labels.append('Alt: %s' % donor)
        if files_ref[d] != 'NA':
            cax.append(cv.cov_from_bam(event_chr, event_pos.min(), event_pos.max(), [files_ref[d]], color_cov='0.0', ax=ax, intron_cnt=True, color_intron_edge='0.0', log=('--log' in sys.argv), intron_filter=intron_filter, return_legend_handle=True, label='Ref: %s' % ref_donors[d])) 
            labels.append('Ref: %s' % ref_donors[d])
        plt.legend(cax, labels, fontsize=10)
        #for vp in var_pos:
        #    ax.plot(vp, 0, 'bo', markersize=1)
        #ylim = ax.get_ylim()
        #xlim = ax.get_xlim()
        #xspan = xlim[1] - xlim[0]
        #yspan = ylim[1] - ylim[0]
        #for vp in var_pos:
        #    ax.arrow(vp, yspan * 0.2, 0, yspan * -0.15, head_width=0.01*xspan, head_length=0.01*yspan, fc='b') 
        x_range = ax.get_xlim()
        ax.set_ylim([0, ax.get_ylim()[1]])
        axs.clean_axis(ax)
        axs.set_ticks_outer(ax)
   
    ### plot event layout
    ax = fig.add_subplot(gs[rows - 1, 0])
    gl.multiple([exons1, exons2], x_range=x_range, ax=ax)
    axs.clean_axis(ax, left=True, bottom=True)

    plt.tight_layout()

    outbase = os.path.join(PLOTDIR, 'result_exonization_exon_skip_%s_%s_%s' % (results[i, 1], results[i, 0], un.get_ID(results[i, 9], lookup=lookup)))
    print('saving to %s' % outbase)
    plt.savefig('%s%s.pdf' % (outbase, logtag), bbox_inches='tight', format='pdf')
    plt.savefig('%s%s.png' % (outbase, logtag), bbox_inches='tight', format='png')
    plt.close(fig)

