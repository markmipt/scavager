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

redcolor='#FC6264'
bluecolor='#70aed1'
greencolor='#8AA413'

def get_basic_distributions(df):
    mz_array = ((df['calc_neutral_pep_mass'] + df['assumed_charge'] * 1.007276 ) / df['assumed_charge']).values
    rt_exp_array = (df['RT exp']).values
    lengths_array = df['length'].values
    return mz_array, rt_exp_array, lengths_array

def plot_basic_figures(df, df_f, fig):
    mz_array, rt_exp_array, lengths_array = get_basic_distributions(df)
    mz_array_valid, rt_exp_array_valid, lengths_array_valid = get_basic_distributions(df_f)
    fig.add_subplot(1, 3, 1)
    cbins = get_bins(mz_array, mz_array_valid, nbins=100)
    plt.hist(mz_array, bins=cbins, color=redcolor, alpha=0.8, edgecolor='#EEEEEE')
    plt.hist(mz_array_valid, bins=cbins, color=greencolor, alpha=0.8, edgecolor='#EEEEEE')
    plt.ylabel('# of identifications')
    plt.xlabel('PSMs, precursor m/z')
    fig.add_subplot(1, 3, 2)
    cbins = get_bins(rt_exp_array, rt_exp_array_valid, nbins=100)
    plt.hist(rt_exp_array, bins=cbins, color=redcolor, alpha=0.8, edgecolor='#EEEEEE')
    plt.hist(rt_exp_array_valid, bins=cbins, color=greencolor, alpha=0.8, edgecolor='#EEEEEE')
    plt.ylabel('# of identifications')
    plt.xlabel('PSMs, RT experimental')
    fig.add_subplot(1, 3, 3)
    cbins = get_bins(lengths_array, lengths_array_valid, step=1)
    plt.hist(lengths_array, bins=cbins, color=redcolor, alpha=0.8,edgecolor='#EEEEEE')
    plt.hist(lengths_array_valid, bins=cbins, color=greencolor, alpha=0.8,edgecolor='#EEEEEE')
    plt.ylabel('# of identifications')
    plt.xlabel('PSMs, peptide length')

def get_bins(arr1, arr2, nbins=False, step=False):
    tmp = np.concatenate((arr1, arr2))
    minv = tmp.min()
    maxv = tmp.max()
    if nbins:
        return np.linspace(minv, maxv+1, num=nbins)
    elif step:
        return np.arange(minv, maxv+1, step)

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

def plot_outfigures(df, df_f, outfolder, outbasename):
    fig = plt.figure(figsize=(16, 12))
    dpi = fig.get_dpi()
    fig.set_size_inches(2000.0/float(dpi), 2000.0/float(dpi))
    plot_basic_figures(df, df_f, fig)
    plt.grid(color='#EEEEEE')
    plt.tight_layout()
    plt.savefig(path.join(outfolder, outbasename) + '.png')