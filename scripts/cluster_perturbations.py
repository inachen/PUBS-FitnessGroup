#!/bin/python
from __future__ import division
import sys
import os
import argparse
import cPickle as pickle
import numpy as np
import scipy
from scipy import stats
from scipy.cluster.hierarchy import single, linkage, complete, ward, dendrogram
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pylab

NAN_REPLACEMENT = -100 #number to replace np.nan with for downstream analysis, should be absurdly large/small

parser = argparse.ArgumentParser()
parser.add_argument('array_pickles', type=str, nargs='+', help='pickle files encoding numpy arrays of relative fitness values or interactions')
parser.add_argument('--cluster_method', type=str, choices=['all', 'pos'], default='all', help="""mode to use for clustering arrays.
                                         "all" computes pairwise distances between all values of matrices and uses to cluster.
                                         "pos" projects array values of all mutants at a position to a single mean value,
                                         then clusters resulting vectors""")
parser.add_argument('--pert_names', type=str, nargs='+', default=None, help='list of names of perturbations (used for plotting hierarchy)')
parser.add_argument('--data_type', type=str, choices=['rel_fitness','interaction'], default='rel_fitness', help='type of data in input arrays, used for color scheme [choices: rel_fitness or interaction]')
parser.add_argument('--out_plot', type=str, default=None, help='filename for output plot')
args = parser.parse_args()

data_array_list = []

for data_file in args.array_pickles:
    data_array = pickle.load(open(data_file, 'rb'))
    data_array[np.where(data_array==NAN_REPLACEMENT)] = np.nan
    data_array_list.append(data_array)

data_arrays = np.array(data_array_list)
mean_data_array_list = [np.nanmean(x, axis=0) for x in data_array_list]
mean_data_arrays = np.array(mean_data_array_list)


if args.cluster_method is 'all': #calculate distances based on every position that is real in all arrays
    real_indices = np.where(~np.isnan(np.sum(data_arrays, axis=0)))
    real_data_array_list = [x[real_indices] for x in data_array_list]

    distance_matrix = np.empty(shape=(len(real_data_array_list),len(real_data_array_list)))
    for i,x in enumerate(real_data_array_list):
        for j,y in enumerate(real_data_array_list):
            distance_matrix[i][j] = np.linalg.norm(y-x)
else: #calculate distance between position means in all arrays
    real_indices = np.where(~np.isnan(np.sum(mean_data_arrays, axis=0)))
    real_data_array_list = [x[real_indices] for x in mean_data_array_list]

    distance_matrix = np.empty(shape=(len(real_data_array_list),len(real_data_array_list)))
    for i,x in enumerate(real_data_array_list):
        for j,y in enumerate(real_data_array_list):
            distance_matrix[i][j] = np.linalg.norm(y-x)

#determine colormap based on input data
if args.data_type == 'rel_fitness':
    cmap = pylab.cm.YlGnBu_r
    mean_stop = np.nanmin(np.array([np.nanmean(x[0]) for x in data_arrays]))
    vmin = max(-1, mean_stop) #null or lowest
    vmax = 0 #WT
elif args.data_type == 'interaction':
    print "test"
    cmap = pylab.cm.seismic
    max_abs = max(abs(np.nanmin(mean_data_arrays)), abs(np.nanmax(mean_data_arrays)))
    vmin = -max_abs
    vmax = max_abs

clusters = ward(distance_matrix)

# Generate random features and distance matrix.
# Compute and plot first dendrogram.
fig = pylab.figure(figsize=(16,3.5))

#plot dendrogram
ax1 = fig.add_axes([0.04,0.18,0.05,0.50])
dendro = dendrogram(clusters, orientation='right')
ax1.set_xticks([])
ax1.set_yticks([])

# Plot distance matrix.
axmatrix = fig.add_axes([0.09,0.18,0.8,0.50])
idx1 = dendro['leaves']
sorted_distance_matrix = distance_matrix[idx1,:]
sorted_mean_data_array = mean_data_arrays[idx1,:]
if args.pert_names is not None:
    pert_names = np.array(args.pert_names)
    sorted_labels = pert_names[idx1]
else:
    sorted_labels=[]
im = axmatrix.matshow(sorted_mean_data_array, aspect='auto', origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
axmatrix.set_xticks(range(0,76))
axmatrix.xaxis.set_ticks_position('bottom')
axmatrix.set_xticklabels(range(2,78), fontsize=9, y=0.02)
axmatrix.set_yticks(range(len(sorted_labels)))
axmatrix.set_yticklabels(sorted_labels, fontsize=12)
axmatrix.yaxis.tick_right()
axmatrix.yaxis.set_label_position('right')
axmatrix.grid(False)
axmatrix.set_title("Mean Fitness for Perturbations", y=1.35, fontsize=18)
axmatrix.set_ylabel("Perturbations", labelpad=10, fontsize=14)
axmatrix.set_xlabel("Position",labelpad=10, fontsize=14)


# Plot colorbar.
axcolor = fig.add_axes([0.09,0.7,0.8,0.05])
pylab.colorbar(im, cax=axcolor, orientation='horizontal')
axcolor.xaxis.tick_top()
axcolor.xaxis.set_label_position('top')


#fig.set_tight_layout(True)
if args.out_plot is not None:
    fig.savefig(args.out_plot)
else:
    pylab.show()