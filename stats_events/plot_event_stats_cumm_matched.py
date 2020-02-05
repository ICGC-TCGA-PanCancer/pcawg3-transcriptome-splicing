import sys
import os
import scipy as sp
import pickle
import h5py

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.linewidth'] = 2
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.append('/cluster/home/akahles/git/tools/python/viz/')
import axes as axs
sys.path.append('/cluster/home/akahles/git/software/spladder/python')

sys.path.append('..')
from paths import *

BASEDIR_GTEX = os.path.join(BASEDIR_AS, 'alternative_splicing_gtex_matched')
BASEDIR_ICGC = os.path.join(BASEDIR_AS, 'alternative_splicing_icgc_matched')
BASEDIR_ANNO = os.path.join(BASEDIR_AS, 'alternative_splicing_anno')
PLOTDIR      = os.path.join(BASEDIR_AS, 'alternative_splicing/plots/event_stats')
CONF = 2

event_types = ['exon_skip',
               'intron_retention',
               'alt_3prime',
               'alt_5prime',
               'mutex_exons',
               'mult_exon_skip']
if CONF == 3:
    event_types = event_types[:-1]
event_dict = {'exon_skip':'Exon Skip',
              'intron_retention':'Intron Retention',
              'alt_3prime':'Alternative 3\' Site',
              'alt_5prime':'Alternative 5\' Site',
              'mutex_exons':'Mutually Exclusive Exons',
              'mult_exon_skip':'Mult. Exon Skip'}

def main():
    figs = dict()
    figs['stats'] = plt.figure(figsize=(12, 8))
    figs['stats_log'] = plt.figure(figsize=(12, 8))
    figs['stats_full'] = plt.figure(figsize=(12, 8))
    figs['stats_full_log'] = plt.figure(figsize=(12, 8))
    gss = dict()
    gss['stats'] = gridspec.GridSpec(2, 3) #, wspace=0.0, hspace=0.0)
    gss['stats_log'] = gridspec.GridSpec(2, 3) #, wspace=0.0, hspace=0.0)
    gss['stats_full'] = gridspec.GridSpec(2, 3) #, wspace=0.0, hspace=0.0)
    gss['stats_full_log'] = gridspec.GridSpec(2, 3) #, wspace=0.0, hspace=0.0)

    for e, event_type in enumerate(event_types):

        print('Handling %s' % event_type, file=sys.stderr)
        
        ### load events detected in annotation only
        anno = pickle.load(open(os.path.join(BASEDIR_ANNO, 'merge_graphs_%s_C%i.pickle' % (event_type, CONF)), 'r'))
        if isinstance(anno, tuple):
            anno = anno[0]

        ### load annotation index
        is_anno_gtex = pickle.load(open(os.path.join(BASEDIR_GTEX, 'merge_graphs_%s_C%i.anno_only.pickle' % (event_type, CONF)), 'r'))
        is_anno_icgc = pickle.load(open(os.path.join(BASEDIR_ICGC, 'merge_graphs_%s_C%i.anno_only.pickle' % (event_type, CONF)), 'r'))

        ### load confident events
        IN = h5py.File(os.path.join(BASEDIR_GTEX, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, CONF)), 'r')
        idx_conf_gtex = IN['conf_idx'][:]
        IN.close()
        IN = h5py.File(os.path.join(BASEDIR_ICGC, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, CONF)), 'r')
        idx_conf_icgc = IN['conf_idx'][:]
        IN.close()

        ### load filtered events
        #IN = h5py.File(os.path.join(BASEDIR_GTEX, 'merge_graphs_%s_C%i.counts.r10_s50_V10.hdf5' % (event_type, CONF)), 'r')
        #idx_filt_gtex = IN['filter_idx'][:]
        #IN.close()
        #IN = h5py.File(os.path.join(BASEDIR_ICGC, 'merge_graphs_%s_C%i.counts.r10_s50_V10.hdf5' % (event_type, CONF)), 'r')
        #idx_filt_icgc = IN['filter_idx'][:]
        #IN.close()

        ### load psi filtered events
        idx_psi_gtex = pickle.load(open(os.path.join(BASEDIR_GTEX, 'merge_graphs_%s_C%i.counts.hdf5.psi_filt.pickle' % (event_type, CONF)), 'r'))[1]
        idx_psi_icgc = pickle.load(open(os.path.join(BASEDIR_ICGC, 'merge_graphs_%s_C%i.counts.hdf5.psi_filt.pickle' % (event_type, CONF)), 'r'))[1]

        ### plot stats for normal counts (FULL)
        ax = figs['stats_full'].add_subplot(gss['stats_full'][e / 3, e % 3])
        xlabels_full = ['detected', 'confident']
        xlabels_part = ['confident']
        xlabels_full.extend(['dpsi > %.1f' % _ for _ in sorted(idx_psi_gtex.keys())])
        xlabels_part.extend(['dpsi > %.1f' % _ for _ in sorted(idx_psi_gtex.keys())])
        # all confirmed events, further filtered by PSI - GTEX
        data1_gtex = [is_anno_gtex.shape[0], idx_conf_gtex.shape[0]]
        data1_gtex.extend([sp.intersect1d(idx_conf_gtex, idx_psi_gtex[_]).shape[0] for _ in sorted(idx_psi_gtex.keys())])
        data1_gtex = sp.array(data1_gtex)
        lg, = ax.plot(sp.arange(data1_gtex.shape[0]), data1_gtex, '-b', label='GTEx')
        # all annotated confirmed events, further filtered by PSI - GTEX
        data2_gtex = [sp.sum(is_anno_gtex), sp.sum(is_anno_gtex[idx_conf_gtex])]
        data2_gtex.extend([sp.sum(is_anno_gtex[sp.intersect1d(idx_conf_gtex, idx_psi_gtex[_])]) for _ in sorted(idx_psi_gtex.keys())])
        data2_gtex = sp.array(data2_gtex)
        lga, = ax.plot(sp.arange(data2_gtex.shape[0]), data2_gtex, '--b', label='GTEx (anno)')
        # all confirmed events, further filtered by PSI - ICGC
        data1_icgc = [is_anno_icgc.shape[0], idx_conf_icgc.shape[0]]
        data1_icgc.extend([sp.intersect1d(idx_conf_icgc, idx_psi_icgc[_]).shape[0] for _ in sorted(idx_psi_icgc.keys())])
        data1_icgc = sp.array(data1_icgc)
        li, = ax.plot(sp.arange(data1_icgc.shape[0]), data1_icgc, '-r', label='ICGC')
        # all annotated confirmed events, further filtered by PSI - ICGC
        data2_icgc = [sp.sum(is_anno_icgc), sp.sum(is_anno_icgc[idx_conf_icgc])]
        data2_icgc.extend([sp.sum(is_anno_icgc[sp.intersect1d(idx_conf_icgc, idx_psi_icgc[_])]) for _ in sorted(idx_psi_icgc.keys())])
        data2_icgc = sp.array(data2_icgc)
        lia, = ax.plot(sp.arange(data2_icgc.shape[0]), data2_icgc, '--r', label='ICGC (anno)')
        axs.set_ticks_outer(ax)
        axs.clean_axis(ax)
        if e == len(event_types) - 1:
            ax.legend(handles=[lg, lga, li, lia], loc='upper right', frameon=False, fontsize=10)
        ax.set_xticks(list(range(len(xlabels_full))))
        if e < len(event_types) - 3:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(xlabels_full, rotation=90, fontsize=10)
        ax.set_title(event_dict[event_type])
        ax.xaxis.grid(True)


        ### plots stats for log10 counts (FULL)
        ax = figs['stats_full_log'].add_subplot(gss['stats_full_log'][e / 3, e % 3])
        lg, = ax.plot(sp.arange(data1_gtex.shape[0]), sp.log10(data1_gtex + 1), '-b', label='GTEx')
        lga, = ax.plot(sp.arange(data2_gtex.shape[0]), sp.log10(data2_gtex + 1), '--b', label='GTEx (anno)')
        li, = ax.plot(sp.arange(data1_icgc.shape[0]), sp.log10(data1_icgc + 1), '-r', label='ICGC')
        lia, = ax.plot(sp.arange(data2_icgc.shape[0]), sp.log10(data2_icgc + 1), '--r', label='ICGC (anno)')
        axs.set_ticks_outer(ax)
        axs.clean_axis(ax)
        if e == len(event_types) - 1:
            ax.legend(handles=[lg, lga, li, lia], loc='lower left', frameon=False, fontsize=10)
        ax.set_xticks(list(range(len(xlabels_full))))
        if e < len(event_types) - 3:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(xlabels_full, rotation=90, fontsize=10)
        ax.set_title(event_dict[event_type])
        ax.xaxis.grid(True)

        ### plot stats for normal counts (only conf)
        ax = figs['stats'].add_subplot(gss['stats'][e / 3, e % 3])
        lg, = ax.plot(sp.arange(data1_gtex.shape[0] - 1), data1_gtex[1:], '-b', label='GTEx')
        lga, = ax.plot(sp.arange(data2_gtex.shape[0] - 1), data2_gtex[1:], '--b', label='GTEx (anno)')
        ax.plot(sp.arange(data1_icgc.shape[0] - 1), data1_icgc[1:], '-r', label='ICGC')
        ax.plot(sp.arange(data2_icgc.shape[0] - 1), data2_icgc[1:], '--r', label='ICGC (anno)')
        axs.set_ticks_outer(ax)
        axs.clean_axis(ax)
        if e == len(event_types) - 1:
            ax.legend(handles=[lg, lga, li, lia], loc='upper right', frameon=False, fontsize=10)
        ax.set_xticks(list(range(len(xlabels_part))))
        if e < len(event_types) - 3:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(xlabels_part, rotation=90, fontsize=10)
        ax.set_title(event_dict[event_type])
        ax.xaxis.grid(True)

        ### plots stats for log10 counts (only cony)
        ax = figs['stats_log'].add_subplot(gss['stats_log'][e / 3, e % 3])
        lg, = ax.plot(sp.arange(data1_gtex.shape[0] - 1), sp.log10(data1_gtex[1:] + 1), '-b', label='GTEx')
        lga, = ax.plot(sp.arange(data2_gtex.shape[0] - 1), sp.log10(data2_gtex[1:] + 1), '--b', label='GTEx (anno)')
        li, = ax.plot(sp.arange(data1_icgc.shape[0] - 1), sp.log10(data1_icgc[1:] + 1), '-r', label='ICGC')
        lia, = ax.plot(sp.arange(data2_icgc.shape[0] - 1), sp.log10(data2_icgc[1:] + 1), '--r', label='ICGC (anno)')
        axs.set_ticks_outer(ax)
        axs.clean_axis(ax)
        if e == len(event_types) - 1:
            ax.legend(handles=[lg, lga, li, lia], loc='lower left', frameon=False, fontsize=10)
        ax.set_xticks(list(range(len(xlabels_part))))
        if e < len(event_types) - 3:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(xlabels_part, rotation=90, fontsize=10)
        ax.set_title(event_dict[event_type])
        ax.xaxis.grid(True)


    for p in figs:
        figs[p].tight_layout()
        figs[p].savefig(os.path.join(PLOTDIR, 'event_overview_cumm_matched_C%i_%s.pdf' % (CONF, p)), format='pdf', bbox_inches='tight')
        figs[p].savefig(os.path.join(PLOTDIR, 'event_overview_cumm_matched_C%i_%s.png' % (CONF, p)), format='png', bbox_inches='tight')
        plt.close(figs[p])

if __name__ == "__main__":
    main()
