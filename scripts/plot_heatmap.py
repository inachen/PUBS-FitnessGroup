#!/usr/bin/python
from __future__ import division
import os
import sys
import cPickle as pkl
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse

'''Input numpy arrays with amino acid fitness or interaction data, and plot as heatmaps.'''

def plot_fitness_heatmap(data_array, ax, start_pos, minv, maxv, cmap, ylabels=[], 
                         plot_title='', color_ticks=None, color_tick_labels=None):
    '''Given a data array and a plot axis, creates heatmap and color bar and adds to axis'''

    pos = np.arange(start_pos,data_array.shape[1]+start_pos)
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(pos.shape[0]), minor=True)
    ax.set_xticklabels([x if i%5==0 else '' for i,x in enumerate(pos)], minor=True) #show every 5 positions on plots
    ax.set_xlim([-0.5,len(pos)-0.5])
    ax.set_xlabel("Position", labelpad=8, fontsize=16)

    ax.yaxis.set_ticks_position('left') # this one is optional but I still recommend it...
    ax.set_yticks(np.arange(len(ylabels)))
    ax.set_yticklabels(ylabels)
    ax.set_ylabel("Amino Acid", labelpad=0, fontsize=16)

    ax.set_title(plot_title, fontsize=20, y=1.05)
    ax.set_aspect('equal')
    cmap.set_bad('1') #NaNs will be shown in white
    im = ax.imshow(data_array*~np.isnan(data_array), cmap=cmap, interpolation="nearest", vmin=minv, vmax=maxv)

    ax.grid(False)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1.5%", pad=0.1)
    if color_ticks is not None and color_tick_labels is not None:
        cbar = plt.colorbar(im, cax=cax, ticks=color_ticks, orientation='vertical')
        cax.set_yticklabels(color_tick_labels)
    else:
        cbar = plt.colorbar(im, cax=cax, orientation='vertical')
    cbar.ax.tick_params(labelsize=14)

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)
    return newcmap

#run
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Given a pickle containing fitness values or interactions, plots heatmap")
    #necessary input_files
    parser.add_argument('data', type=str, default=None, nargs='+', help='pickle file(s) containing numpy array with amino fitness or interaction data')
    parser.add_argument('--aa_index', type=str, default=None, help='pickle file encoding a dictionary of amino acids and corresponding indices used to label plot')
    #running options
    parser.add_argument('--start_pos', type=int, default=2, help='first position on the plot')
    parser.add_argument('--plot_titles', type=str, default=None, nargs='+', help='Label of sample, used when making multiple plots [e.g. "DMSO"]')
    parser.add_argument('--data_type', type=str, default='fitness', choices=('fitness', 'interaction', 'other'), help='data type encoded by the pickle, used for setting color maps [default: fitness (YlGnBu colormap)]')
    parser.add_argument('--min', type=float, default=None, help='lowest value on the colormap')
    parser.add_argument('--center', type=float, default=None, help='central value on the colormap')
    parser.add_argument('--max', type=float, default=None, help='largest value on the colormap')
    parser.add_argument('-o', '--out_plot', type=str, default=None, help='output file to save plot')
    parser.add_argument('--nan_standins', type=int, default=-100, help='integer NaNs are cast to [default:-100]')
    args = parser.parse_args()

    y_axis_labels=[]
    if args.aa_index is not None:
        aa_index = pkl.load(open(args.aa_index, 'rb'))
        y_axis_labels = sorted(aa_index.keys(), key=lambda x: aa_index[x])

    #set color bar parameters
    minv=None
    maxv=None
    color_ticks = None
    color_tick_labels = None
    cmap=plt.cm.Spectral_r

    if args.data_type == 'fitness':
        minv=-1
        maxv=0
        color_ticks = [-1, 0]
        color_tick_labels = ['Null', 'WT']
        cmap=plt.cm.YlGnBu_r
    elif args.data_type == 'interaction':
        minv=-1
        maxv=1
        color_ticks = [-1, 0, 1]
        color_tick_labels = ['Negative', 'Neutral', 'Positive']
        cmap=plt.cm.seismic_r

    if args.min is not None:
        minv=args.min
        color_ticks=None
        color_tick_labels=None
    if args.max is not None:
        maxv=args.max
        color_ticks=None
        color_tick_labels=None


    fig = plt.figure(figsize=(18,5*len(args.data)))
    fig.subplots_adjust(hspace=.3)

    for i,data_file in enumerate(args.data):
        data_array = pkl.load(open(data_file, 'rb'))
        data_array[np.where(data_array == args.nan_standins)]=np.nan #replace nan standins with nans
        
        #set special color bar parameters for subplot 
        if minv is None:
            minv = np.nanmin(data_array)
        if maxv is None:
            maxv = np.nanmax(data_array)
        if args.center is not None:
            new_cmap=shiftedColorMap(cmap, start=0, midpoint=1-maxv/(maxv+abs(minv)), stop=1)
        else:
            new_cmap=cmap

        #set plot title
        plot_title=''
        if args.plot_titles is not None and i < len(args.plot_titles):
            plot_title = args.plot_titles[i]

        #add subplot
        ax = fig.add_subplot(len(args.data),1,i+1)
        plot_fitness_heatmap(data_array, ax, start_pos=args.start_pos, minv=minv, maxv=maxv, 
                             cmap=new_cmap, ylabels=y_axis_labels, plot_title=plot_title,
                             color_ticks=color_ticks, color_tick_labels=color_tick_labels)

    plt.tight_layout()

    #save or show plot
    if args.out_plot is not None:
        plt.savefig(args.out_plot, dpi=300)
    else:
        plt.show()
