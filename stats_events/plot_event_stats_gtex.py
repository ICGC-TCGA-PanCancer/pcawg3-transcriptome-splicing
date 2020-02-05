import sys
import os
import scipy as sp
import pickle
import h5py

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.append('/cluster/home/akahles/git/tools/python/viz/')
import axes as axs
sys.path.append('/cluster/home/akahles/git/software/spladder/python')

sys.path.append('..')
from paths import *

BASEDIR_GTEX = os.path.join(BASEDIR_AS, 'alternative_splicing')
BASEDIR_ANNO = os.path.join(BASEDIR_AS, 'alternative_splicing_anno')
PLOTDIR      = os.path.join(BASEDIR_AS, 'alternative_splicing/plots/event_stats')
CONF = 3

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

def plot_full(ax, vals, ylabel, title):
    ### detected / detected conf / detected high-conf
    cax = ax.bar([2, 5, 8], vals[1], color=['b', 'g', 'm'], alpha=0.5, width=0.8, linewidth=0)
    axs.label_bars(ax, cax, vals[1], rotation=90)
    ### annotated fraction of the counts above
    cax = ax.bar([3, 6, 9], vals[2], color=['b', 'g', 'm'], alpha=0.9, width=0.8, linewidth=0)
    axs.label_bars(ax, cax, vals[2], rotation=90)
    ### high-conf novel events
    cax = ax.bar([11], vals[3], color=['c'], alpha=0.9, width=0.8, linewidth=0)
    axs.label_bars(ax, cax, vals[3], rotation=90)
    ### annotated
    cax = ax.bar([0], vals[0], color='r', alpha=0.5, width=0.8, linewidth=0)
    axs.label_bars(ax, cax, vals[0], rotation=90)

    ax.set_ylabel(ylabel)
    ax.set_xlabel('')
    ax.set_xticks(sp.array([0, 2, 3, 5, 6, 8, 9, 11], dtype='float') + 0.4)
    ax.set_xticklabels(['annotated', 'detected', 'detected\n(anno)', 'detected\n(conf)', 'detected\n(conf+anno)', 'detected\n(h_conf)', 'detected\n(h_conf+anno)', 'novel\n(h_conf)'], rotation=90)
    ax.set_xlim([-0.2, ax.get_xlim()[1]])
    ax.set_title(title)

    axs.set_ticks_outer(ax)
    axs.clean_axis(ax)

def plot_part(ax, vals, ylabel, title):
    ### detected / detected conf / detected high-conf
    cax = ax.bar([2, 5], vals[1], color=['g', 'm'], alpha=0.5, width=0.8, linewidth=0)
    axs.label_bars(ax, cax, vals[1], rotation=90)
    ### annotated fraction of the counts above
    cax = ax.bar([3, 6], vals[2], color=['g', 'm'], alpha=0.9, width=0.8, linewidth=0)
    axs.label_bars(ax, cax, vals[2], rotation=90)
    ### high-conf novel events
    cax = ax.bar([8], vals[3], color=['c'], alpha=0.9, width=0.8, linewidth=0)
    axs.label_bars(ax, cax, vals[3], rotation=90)
    ### annotated
    cax = ax.bar([0], vals[0], color='r', alpha=0.5, width=0.8, linewidth=0)
    axs.label_bars(ax, cax, vals[0], rotation=90)

    ax.set_ylabel(ylabel)
    ax.set_xlabel('')
    ax.set_xticks(sp.array([0, 2, 3, 5, 6, 8], dtype='float') + 0.4)
    ax.set_xticklabels(['annotated', 'detected\n(conf)', 'detected\n(conf+anno)', 'detected\n(h_conf)', 'detected\n(h_conf+anno)', 'novel\n(h_conf)'], rotation=90)
    ax.set_xlim([-0.2, ax.get_xlim()[1]])
    ax.set_title(title)

    axs.set_ticks_outer(ax)
    axs.clean_axis(ax)



def main():
    figs = dict()
    figs['stats'] = plt.figure(figsize=(10, 12))
    figs['stats_log'] = plt.figure(figsize=(10, 12))
    figs['stats_full'] = plt.figure(figsize=(10, 12))
    figs['stats_full_log'] = plt.figure(figsize=(10, 12))
    gss = dict()
    gss['stats'] = gridspec.GridSpec(3, 2) #, wspace=0.0, hspace=0.0)
    gss['stats_log'] = gridspec.GridSpec(3, 2) #, wspace=0.0, hspace=0.0)
    gss['stats_full'] = gridspec.GridSpec(3, 2) #, wspace=0.0, hspace=0.0)
    gss['stats_full_log'] = gridspec.GridSpec(3, 2) #, wspace=0.0, hspace=0.0)

    for e, event_type in enumerate(event_types):

        print('Handling %s' % event_type, file=sys.stderr)
        
        ### load events detected in annotation only
        anno = pickle.load(open(os.path.join(BASEDIR_ANNO, 'merge_graphs_%s_C%i.pickle' % (event_type, CONF)), 'r'))
        if isinstance(anno, tuple):
            anno = anno[0]

        ### load annotation index
        is_anno = pickle.load(open(os.path.join(BASEDIR, 'merge_graphs_%s_C%i.anno_only.pickle' % (event_type, CONF)), 'r'))

        ### load confident events
        IN = h5py.File(os.path.join(BASEDIR, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, CONF)), 'r')
        idx_conf = IN['conf_idx'][:]
        IN.close()

        ### load filtered events
        IN = h5py.File(os.path.join(BASEDIR, 'merge_graphs_%s_C%i.counts.r10_s50_V10.hdf5' % (event_type, CONF)), 'r')
        idx_filt = IN['filter_idx'][:]
        IN.close()


        ### plot stats for normal counts (FULL)
        ax = figs['stats_full'].add_subplot(gss['stats_full'][e / 2, e % 2])
        vals = dict()
        vals[0] = [anno.shape[0]]
        vals[1] = [is_anno.shape[0], idx_conf.shape[0], idx_filt.shape[0]]
        vals[2] = [sp.sum(is_anno), sp.sum(is_anno[idx_conf]), sp.sum(is_anno[idx_filt])]
        vals[3] = [sp.sum(~is_anno[idx_filt])]

        plot_full(ax, vals, '# of events', event_dict[event_type])

        ### plots stats for log10 counts (FULL)
        ax = figs['stats_full_log'].add_subplot(gss['stats_full_log'][e / 2, e % 2])
        vals = dict()
        vals[0] = sp.log10([anno.shape[0]])
        vals[1] = sp.log10([is_anno.shape[0], idx_conf.shape[0], idx_filt.shape[0]])
        vals[2] = sp.log10([sp.sum(is_anno), sp.sum(is_anno[idx_conf]), sp.sum(is_anno[idx_filt])])
        vals[3] = sp.log10([sp.sum(~is_anno[idx_filt])])

        plot_full(ax, vals, '# of events (log 10)', event_dict[event_type])

        ### plot stats for normal counts (only conf)
        ax = figs['stats'].add_subplot(gss['stats'][e / 2, e % 2])
        vals = dict()
        vals[0] = [anno.shape[0]]
        vals[1] = [idx_conf.shape[0], idx_filt.shape[0]]
        vals[2] = [sp.sum(is_anno[idx_conf]), sp.sum(is_anno[idx_filt])]
        vals[3] = [sp.sum(~is_anno[idx_filt])]

        plot_part(ax, vals, '# of events', event_dict[event_type])

        ### plots stats for log10 counts (only cony)
        ax = figs['stats_log'].add_subplot(gss['stats_log'][e / 2, e % 2])
        vals = dict()
        vals[0] = sp.log10([anno.shape[0]])
        vals[1] = sp.log10([idx_conf.shape[0], idx_filt.shape[0]])
        vals[2] = sp.log10([sp.sum(is_anno[idx_conf]), sp.sum(is_anno[idx_filt])])
        vals[3] = sp.log10([sp.sum(~is_anno[idx_filt])])

        plot_part(ax, vals, '# of events (log 10)', event_dict[event_type])

    for p in figs:
        figs[p].tight_layout()
        figs[p].savefig(os.path.join(PLOTDIR, 'event_overview_gtex_C%i_%s.pdf' % (CONF, p)), format='pdf', bbox_inches='tight')
        figs[p].savefig(os.path.join(PLOTDIR, 'event_overview_gtex_C%i_%s.png' % (CONF, p)), format='png', bbox_inches='tight')
        plt.close(figs[p])

if __name__ == "__main__":
    main()
