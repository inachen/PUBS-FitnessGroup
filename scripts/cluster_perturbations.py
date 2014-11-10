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

#Functions



parser = argparse.ArgumentParser()
parser.add_argument('rel_aa_fitnesses', type=str, nargs='+', help='pickle files encoding numpy arrays of relative fitness values')
parser.add_argument('--cluster_method', type=str, choices=['all', 'pos'], help="""mode to use for clustering fitness matrices.
                                         "all" computes pairwise distances between all values of matrices and uses to cluster.
                                         "pos" projects fitness values of all mutants at a position to a single mean fitness value,
                                         then clusters resulting fitness vectors""")
parser.add_argument('--pert_names', type=str, nargs='+', help='list of names of perturbations (used for plotting hierarchy)')
args = parser.parse_args()

fitness_array_list = []

for fitness_file in args.rel_aa_fitnesses:
    fitness_array = pickle.load(open(fitness_file, 'rb'))
    fitness_array_list.append(fitness_array)

fitness_arrays = np.array(fitness_array_list)
fitness_arrays[np.where(fitness_arrays==NAN_REPLACEMENT)] = np.nan


real_indices = np.where(~np.isnan(np.sum(fitness_arrays, axis=0)))
real_fitness_array_list = [x[real_indices] for x in fitness_array_list]

distance_matrix = np.empty(shape=(len(real_fitness_array_list),len(real_fitness_array_list)))
for i,x in enumerate(real_fitness_array_list):
    for j,y in enumerate(real_fitness_array_list):
        distance_matrix[i][j] = np.linalg.norm(y-x)


clusters = ward(distance_matrix)
'''
header = args.perturbation_names

fig = pylab.figure()

# Compute and plot second dendrogram.
ax1 = fig.add_axes([0.3,.725,0.6,0.2])
dendro = dendrogram(clusters)
idx1 = dendro['leaves']
sorted_labels = [header[x] for x in idx1]
#sorted_colors = [dendro['color_list'][x] for x in idx1]
#print sorted_colors
ax1.set_xticklabels(sorted_labels, fontsize=8, rotation='vertical')
ax1.set_yticks([])

# Plot distance matrix.
axmatrix = fig.add_axes([0.3,0.05,0.6,0.6])
sorted_matrix = np.empty_like(rmsd_matrix)
for i,x in enumerate(idx1):
    for j,y in enumerate(reversed(idx1)):
        sorted_rmsd_matrix[i][j] = rmsd_matrix[x][y]
sorted_rmsd_matrix = sorted_rmsd_matrix.T
norm = colors.Normalize(sorted_rmsd_matrix.min(), sorted_rmsd_matrix.max())
im = axmatrix.matshow(sorted_rmsd_matrix, aspect='auto', origin='lower', cmap=pylab.cm.gist_gray, norm=norm)
axmatrix.set_xticks([])
axmatrix.set_yticks([])
axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
pylab.show()

'''
'''
#Read in files
if args.aa_fitness is not None:
    amino_to_number_dict = pickle.load(open(args.aa_index, "rb"))
    aa_fitness_matrices = []

    for aa_fitness_pickle in args.aa_fitness:
        print aa_fitness_pickle
        aa_fitness_dict = pickle.load(open(aa_fitness_pickle, "rb"))
        aa_fitness_matrix = make_aa_fitness_matrix(aa_fitness_dict, amino_to_number_dict)
        aa_fitness_matrices.append(aa_fitness_matrix)

    matrix1 = average_matrices(aa_fitness_matrices[:3])
    matrix2 = average_matrices(aa_fitness_matrices[3:])

    flat_mat1 = np.ndarray.flatten(matrix1)
    flat_mat2 = np.ndarray.flatten(matrix2)
    filtered_flat_mat1 = np.array([x for i,x in enumerate(flat_mat1) if not np.isnan(x) and not np.isnan(flat_mat2[i])])
    filtered_flat_mat2 = np.array([x for i,x in enumerate(flat_mat2) if not np.isnan(x) and not np.isnan(flat_mat1[i])])
    slope, intercept, r_squ, p = regression(filtered_flat_mat1, filtered_flat_mat2)
    print r_squ,slope
    new_y = slope*filtered_flat_mat1 + intercept
    #print new_y
    plt.scatter(np.ndarray.flatten(matrix1),np.ndarray.flatten(matrix2))
    plt.plot(filtered_flat_mat1, new_y, color='r')
    plt.title("Day 1 vs Day 2 Fitness Comparisons", y=1.05)
    plt.xlabel("Day 1 Fitness", labelpad=10)
    plt.ylabel("Day 2 Fitness", labelpad=10)
    plt.show()
    sys.exit()


    if len(aa_fitness_matrices) == 2:
        if args.average:
            aa_fitness_matrix = average_matrices(aa_fitness_matrices)
        elif args.subtract:
            aa_fitness_matrix = subtract_matrices(aa_fitness_matrices[0], aa_fitness_matrices[1])
        else:
            print "No command given for handling multiple matrices. Only first matrix will be considered."
            aa_fitness_matrix = aa_fitness_matrices[0]
    elif len(aa_fitness_matrices) > 2:
        if args.average and not args.subtract:
            aa_fitness_matrix = average_matrices(aa_fitness_matrices)
        elif args.subtract:
            avg_matrix = average_matrices(aa_fitness_matrices[:-1])
            aa_fitness_matrix = subtract_matrices(avg_matrix, aa_fitness_matrices[-1])
    else:
        aa_fitness_matrix = aa_fitness_matrices[0]

    if args.seq_entropy:
        entropy_values = calculate_sequence_entropy(aa_fitness_matrix, amino_to_number_dict)
        print "attribute: entropy\nrecipient: residues"
        for i, entropy in enumerate(entropy_values):
            print "\t:%d.A\t%f" % (i+1, entropy)

    if args.plot_heatmaps:
        column_labels = [x for x,y in sorted(amino_to_number_dict.items(), key=lambda x: x[1])]
        row_labels = sorted(aa_fitness_dict.keys())
        plt.figure(figsize=(20,6), dpi=72)
        if args.seq_entropy:
            ax1 = plt.subplot(211)
        else:
            ax1 = plt.subplot(111)
        if len(args.aa_fitness) > 1 and args.subtract:
            plt.title("Fitness Difference", y=1.05)
            vmax = np.nanmax(aa_fitness_matrix)
            vmin = np.nanmin(aa_fitness_matrix)
            new_cmap = shiftedColorMap(plt.cm.seismic, start=0, midpoint=1 - vmax/(vmax + abs(vmin)), stop=1)
            heatmap = plt.imshow(aa_fitness_matrix, cmap=new_cmap, interpolation="nearest")
        else:
            plt.title("Fitness", y=1.05)
            #vmax = np.nanmax(aa_fitness_matrix)
            #vmin = np.nanmin(aa_fitness_matrix)
            #new_cmap = shiftedColorMap(plt.cm.seismic, start=0, midpoint=1 - vmax/(vmax + abs(vmin)), stop=1)
            heatmap = plt.imshow(aa_fitness_matrix, cmap=plt.cm.YlGnBu_r, vmax=0,vmin=-0.36, interpolation="nearest")
            #heatmap = plt.imshow(aa_fitness_matrix, cmap=new_cmap, interpolation="nearest")
        #heatmap = plt.imshow(aa_fitness_matrix, cmap=plt.cm.seismic,vmax=0.3,vmin=-0.3, interpolation="nearest")
        ax1.set_xticks(np.arange(row_labels[0], row_labels[-1], 5))
        ax1.set_yticks(np.arange(aa_fitness_matrix.shape[0]), minor=False)
        plt.gca().set_xlim((-0.5, len(row_labels) - 0.5))
        plt.gca().set_ylim((-0.5, len(column_labels) - 0.5))
        ax1.invert_yaxis()
        ax1.set_xticklabels(scipy.arange(row_labels[0] + 1, row_labels[-1] + 1, 5))
        ax1.set_yticklabels(column_labels)

        plt.xlabel("Sequence Position", labelpad=10)
        plt.ylabel("Amino Acid", labelpad=10)
        plt.colorbar(orientation='vertical', shrink=.68, pad=0.01)
        plt.grid(False)
        if args.seq_entropy:
            ax2 = plt.subplot(221)
            entropy_scale = plt.imshow(aa_fitness_matrix, cmap=plt.cm.seismic,vmax=0.3,vmin=-0.3, interpolation="nearest")
        plt.tight_layout()
        plt.show()

'''
'''#Plot heatmaps of barcode counts for each amino acid and position
column_labels = [x for x,y in sorted(amino_to_number_dict.items(), key=lambda x: x[1])]
row_labels = sorted(aa_fitness_scores.keys())
ax1 = plt.subplot(211) #this plot will use a linear color scale
heatmap = plt.imshow(aa_fitness_matrix, cmap=plt.cm.YlGnBu, interpolation="nearest")

ax1.set_xticks(np.arange(row_labels[0], row_labels[-1], 5))
#ax1.set_yticks(np.arange(barcode_count_aa_heatmap.shape[0]), minor=False)
plt.gca().set_xlim((-0.5, len(row_labels) - 0.5))
plt.gca().set_ylim((-0.5, len(column_labels) - 0.5))
ax1.invert_yaxis()
ax1.set_xticklabels(scipy.arange(row_labels[0] + 2, row_labels[-1] + 2, 5))
ax1.set_yticklabels(column_labels)

plt.xlabel("Amino Acid Position")
plt.ylabel("Amino Acid")
plt.title("Amino Acid Fitness")
plt.colorbar(orientation='vertical', shrink=0.75, pad=0.01)
plt.tight_layout()
plt.show()
'''
'''
ax2 = plt.subplot(212) #this plot will use a log color scale
heatmap = plt.imshow(barcode_count_aa_heatmap, norm=LogNorm(vmin=1, vmax=131), cmap=plt.cm.YlGnBu, interpolation="nearest")

ax2.set_xticks(np.arange(row_labels[0], row_labels[-1], 5))
ax2.set_yticks(np.arange(barcode_count_aa_heatmap.shape[0]), minor=False)
plt.gca().set_xlim((-0.5, len(row_labels) - 0.5))
plt.gca().set_ylim((-0.5, len(column_labels) - 0.5))
ax2.invert_yaxis()
ax2.set_xticklabels(scipy.arange(row_labels[0] + 2, row_labels[-1] + 2, 5))
ax2.set_yticklabels(column_labels)

plt.xlabel("Amino Acid Position")
plt.ylabel("Amino Acid")
plt.title("Barcode counts per amino acid at each sequence position (log scale)")
plt.colorbar(orientation='vertical', shrink=0.75, pad=0.01)

plt.tight_layout()
plt.show()


#Plot heatmap of deviation from expected proportion of codon at a given aa position
column_labels = [x for x,y in sorted(codon_to_number_dict.items(), key=lambda x: x[1])]
row_labels = sorted_seq_pos

ax = plt.subplot()
heatmap = plt.imshow(barcode_deviation_from_expected, norm=Normalize(vmin=-0.1, vmax=0.1), cmap=plt.cm.seismic, interpolation="nearest")
ax.set_xticks(np.arange(row_labels[0], row_labels[-1], 5))
ax.set_yticks(np.arange(barcode_deviation_from_expected.shape[0]), minor=False)
plt.gca().set_xlim((-0.5, len(row_labels) - 0.5))
plt.gca().set_ylim((-0.5, len(column_labels) - 0.5))

ax.invert_yaxis()
ax.set_xticklabels(scipy.arange(row_labels[0] + 2, row_labels[-1] + 2, 5))
ax.set_yticklabels(column_labels)

plt.tick_params(axis='both', which='major', labelsize=10)

plt.xlabel("Codon Position")
plt.ylabel("Codon")
plt.title("Deviation of proportion barcode at codon-position from expectation with no bias")

plt.colorbar(orientation='vertical', shrink=0.75, pad=0.01)

plt.tight_layout()
plt.show()
'''

'''
#Plot GC count histograms
ax = plt.subplot()
bins = np.arange(0.0,1.1,0.1)
ax.hist(codon_gc_perc, bins=[0, 0.25, 0.75, 1], color='black', label='Codons')
ax.hist(barcode_gc_perc, bins=bins, color='grey', label='Barcodes')
bin_centers = 0.5*(bins[1:]+bins[:-1])
ax.set_xticks(bin_centers)
ax.set_xticklabels(bins)
plt.gca().set_xlim(-0.1,1.1)
plt.title("GC content of barcodes and codons")
plt.xlabel("GC percent")
plt.ylabel("Count")
plt.legend(loc=2)
plt.show()
'''