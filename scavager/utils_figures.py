from __future__ import division
from scipy.stats import scoreatpercentile
from os import path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
try:
    import seaborn
    seaborn.set(rc={'axes.facecolor':'#ffffff'})
    seaborn.set_style('whitegrid')
except ImportError:
    pass

redcolor = '#FC6264'
bluecolor = '#70aed1'
greencolor = '#8AA413'

def get_basic_distributions(df):
    mz_array = ((df['calc_neutral_pep_mass'] + df['assumed_charge'] * 1.007276 ) / df['assumed_charge']).values
    rt_exp_array = (df['RT exp']).values
    lengths_array = df['length'].values
    return mz_array, rt_exp_array, lengths_array

def get_descriptor_array(df, df_f, dname):
    array_t = df[~df['decoy']][dname].values
    array_d = df[df['decoy']][dname].values
    array_v = df_f[dname].values
    return array_t, array_d, array_v

def plot_hist_basic(array_all, array_valid, fig, subplot_max_x, subplot_i, xlabel, ylabel='# of identifications', bin_size_one=False):
    fig.add_subplot(subplot_max_x, 3, subplot_i)
    cbins = get_bins((array_all, array_valid), bin_size_one)
    plt.hist(array_all, bins=cbins, color=redcolor, alpha=0.8, edgecolor='#EEEEEE')
    plt.hist(array_valid, bins=cbins, color=greencolor, alpha=0.8, edgecolor='#EEEEEE')
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)

def plot_basic_figures(df, df_f, fig, subplot_max_x, subplot_start, idtype):
    mz_array, rt_exp_array, lengths_array = get_basic_distributions(df)
    mz_array_valid, rt_exp_array_valid, lengths_array_valid = get_basic_distributions(df_f)

    plot_hist_basic(mz_array, mz_array_valid, fig, subplot_max_x, subplot_i=subplot_start, \
                     xlabel='%s, precursor m/z' % (idtype, ))
    subplot_start += 1
    plot_hist_basic(rt_exp_array, rt_exp_array_valid, fig, subplot_max_x, subplot_i=subplot_start, \
                     xlabel='%s, RT experimental' % (idtype, ))
    subplot_start += 1
    plot_hist_basic(lengths_array, lengths_array_valid, fig, subplot_max_x, subplot_i=subplot_start, \
                     xlabel='%s, peptide length' % (idtype, ), bin_size_one=True)

def plot_protein_figures(df, df_f, fig, subplot_max_x, subplot_start):
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='LOG10_NSAF'), fig, subplot_max_x, subplot_start, xlabel='LOG10(NSAF)')
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='sq'), fig, subplot_max_x, subplot_start+1, xlabel='sequence coverage')

def plot_hist_descriptor(inarrays, fig, subplot_max_x, subplot_i, xlabel, ylabel='# of identifications'):
    array_t, array_d, array_v = inarrays
    ax = fig.add_subplot(subplot_max_x, 3, subplot_i)
    cbins, width = get_bins_for_descriptors((array_t, array_d, array_v))
    H1, _ = np.histogram(array_d, bins=cbins)
    H2, _ = np.histogram(array_t, bins=cbins)
    H3, _ = np.histogram(array_v, bins=cbins)
    ax.bar(cbins[:-1], H1, width, align='center',color=redcolor, alpha=0.4, edgecolor='#EEEEEE')
    ax.bar(cbins[:-1], H2, width, align='center',color=bluecolor, alpha=0.4, edgecolor='#EEEEEE')
    ax.bar(cbins[:-1], H3, width, align='center',color=greencolor, alpha=1, edgecolor='#EEEEEE')
    cbins=np.append(cbins[0],cbins)
    H1 = np.append(np.append(0,H1),0)
    H2 = np.append(np.append(0,H2),0)
    cbins -= width / 2
    ax.step(cbins, H2, where='post', color=bluecolor, alpha=0.8)
    ax.step(cbins, H1, where='post', color=redcolor, alpha=0.8)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    if 'mass shift' in xlabel:
        plt.xlim(-1.5, 1.5)
        ax.set_xticks([-1, 0, 1])
        ax.set_xticklabels(['NA', 'unmodified', 'modified'])
    elif width == 1.0:
        ax.set_xticks(np.arange(int(cbins[0]), cbins[-1], 1))
        fig.canvas.draw()

def calc_max_x_value(df, df_proteins):
    cnt = 6 # number of basic figures 
    peptide_columns = set(df.columns)
    features_list = ['massdiff_ppm', 'RT diff', 'num_missed_cleavages', 'assumed_charge', 'log_score']
    for feature in features_list:
        if feature in peptide_columns:
            cnt += 1
    for feature in peptide_columns:
        if feature.startswith('mass shift'):
            cnt += 1
    if len(set(df['massdiff_int'])) > 1:
        cnt += 1
    if 'LOG10_NSAF' in df_proteins.columns:
        cnt += 2 # add for NSAF and sequence coverage
    return cnt // 3 + (1 if (cnt % 3) else 0)

    


def plot_descriptors_figures(df, df_f, fig, subplot_max_x, subplot_start):
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='massdiff_ppm'), fig, subplot_max_x, subplot_start, xlabel='precursor mass difference, ppm')
    subplot_start += 1
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='RT diff'), fig, subplot_max_x, subplot_start, xlabel='RT difference, min')
    subplot_start += 1
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='num_missed_cleavages'), fig, subplot_max_x, subplot_start, xlabel='missed cleavages')
    subplot_start += 1
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='assumed_charge'), fig, subplot_max_x, subplot_start, xlabel='precursor charge')
    subplot_start += 1
    if len(set(df['massdiff_int'])) > 1:
        plot_hist_descriptor(get_descriptor_array(df, df_f, dname='massdiff_int'), fig, subplot_max_x, subplot_start, xlabel='isotope mass difference, Da')
        subplot_start += 1
    for df_col in df.columns:
        if df_col.startswith('mass shift'):
            plot_hist_descriptor(get_descriptor_array(df, df_f, dname=df_col), fig, subplot_max_x, subplot_start, xlabel=df_col)
            subplot_start += 1
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='log_score'), fig, subplot_max_x, subplot_start, xlabel='LOG10(ML score)')
    subplot_start += 1
  
def get_bins(inarrays, bin_size_one=False):
    tmp = np.concatenate(inarrays)
    minv = tmp.min()
    maxv = tmp.max()
    if bin_size_one:
        return np.arange(minv, maxv+1, 1)
    else:
        return np.linspace(minv, maxv+1, num=100)

def get_bins_for_descriptors(inarrays):
    tmp = np.concatenate(inarrays)
    minv = tmp.min()
    maxv = tmp.max()
    if len(set(tmp)) <= 15:
        return np.arange(minv, maxv + 2, 1.0), 1.0
    binsize = get_fdbinsize(tmp)
    if binsize < float(maxv - minv) / 300:
        binsize = float(maxv - minv) / 300
    lbin_s = scoreatpercentile(tmp, 1.0)
    lbin = minv
    if lbin_s and abs((lbin - lbin_s) / lbin_s) > 1.0:
        lbin = lbin_s * 1.05
    rbin_s = scoreatpercentile(tmp, 99.0)
    rbin = maxv
    if rbin_s and abs((rbin - rbin_s) / rbin_s) > 1.0:
        rbin = rbin_s * 1.05
    rbin += 1.5 * binsize
    return np.arange(lbin, rbin + binsize, binsize), binsize

def get_fdbinsize(data_list):
    """Calculates the Freedman-Diaconis bin size for
    a data set for use in making a histogram
    Arguments:
    data_list:  1D Data set
    Returns:
    optimal_bin_size:  F-D bin size
    """
    data_list = np.sort(data_list)
    upperquartile = scoreatpercentile(data_list, 75)
    lowerquartile = scoreatpercentile(data_list, 25)
    iqr = upperquartile - lowerquartile
    optimal_bin_size = 2. * iqr / len(data_list) ** (1. / 3.)
    return optimal_bin_size

def plot_outfigures(df, df_f, df_peptides, df_peptides_f, outfolder, outbasename, df_proteins, df_proteins_f):
    fig = plt.figure(figsize=(16, 12))
    dpi = fig.get_dpi()
    fig.set_size_inches(3000.0/float(dpi), 3000.0/float(dpi))
    subplot_max_x = calc_max_x_value(df, df_proteins)
    descriptor_start_index = 7
    plot_basic_figures(df, df_f, fig, subplot_max_x, 1, 'PSMs')
    plot_basic_figures(df_peptides, df_peptides_f, fig, subplot_max_x, 4, 'peptides')
    if 'LOG10_NSAF' in df_proteins.columns:
        plot_protein_figures(df_proteins, df_proteins_f, fig, subplot_max_x, 7)
        descriptor_start_index += 2
    plot_descriptors_figures(df, df_f, fig, subplot_max_x, descriptor_start_index)
    plt.grid(color='#EEEEEE')
    plt.tight_layout()
    plt.savefig(path.join(outfolder, outbasename) + '.png')
    