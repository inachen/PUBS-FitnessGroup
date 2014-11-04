#!/bin/python
import os
import sys
import cPickle as pickle
from scipy.stats import wald
import scipy
import numpy as np
import argparse


def average(l, variances=None):
    '''Calculates simple average of list. If variances for each value are provided, also returns combined variance'''
    mean = 0
    var = 0
    n = len(l)
    if n > 0:
        mean = sum(l)/n
        if variances is not None and len(variances) == n:
            var = sum(variances)/n**2
            return mean, var
        else:
            return mean
    return None

def weighted_average(l, variances):
    '''Calculates the inverse variance weighted average of list and returns average with variance'''
    vals = np.array(l)
    variances = np.array(variances)
    n = len(l)
    if n > 0 and n == len(variances):
        new_var = float(1/np.sum(1/variances**2))
        mean = float(new_var * np.sum(vals/variances**2))
        return new_var*mean
    return None

def rna_to_dna(x):
    '''Converts any RNA sequence to a DNA sequence'''
    return x.replace('U','T').replace('u','t')

def dna_to_rna(x):
    '''Converts any DNA sequence to a RNA sequence'''
    return x.replace('T','U').replace('t','u')

def codon_to_barcodes(pos, codon, allele_dict):
    '''Returns a list of barcodes that map to a specific codon/position combination'''
    return [x for x,y in allele_dict.items() if y[1]==codon and y[0]==pos]

def aa_to_codons(aa, translate_dict):
    '''Returns a list of codons that encode a specific amino acid (or stop)'''
    return [k for k,v in translate_dict.items() if v==aa]

def syn_codons(codon, translate_dict):
    '''Returns a list of codons synonymous to a given codon'''
    return [x for x in aa_to_codons(translate_dict[codon], translate_dict) if x != codon]

def codon_fitness_from_barcodes(barcode_fitness, allele_dict, wt_codon_dict, translate_dict):
    '''Calculates fitness scores for all codons by averaging over fitness values for barcodes and
    returns as a dictionary of positions, then codons, where the corresponding value is just a fitness
    score. Optionally 'remove_outliers=True' can be specified to enable outlier detection'''
    codon_fitness_scores = {}

    for pos in wt_codon_dict.keys():
        codon_fitness_scores.setdefault(pos,{})
        for codon in translate_dict.keys():
            barcode_fitnesses = np.array(sorted([x for x in codon_to_barcodes(pos, codon, allele_dict) if not np.isnan(x)]))
            codon_fitness[pos][codon] = float(np.nanmean(barcode_fitnesses))

    return codon_fitness_scores

def calculate_aa_fitness(codon_fitness_scores, tp_values, wt_counts, total_read_counts, translate_dict, aa_index):
    '''Calculates the fitness of amino acids as the mean fitness of synonymous codons. Returns a dictionary
    with positions as keys, then amino acids (or stop), then the fitness value.'''
    aa_fitness_scores = {}

    for pos in codon_fitness_scores.keys():
        aa_fitness_scores.setdefault(pos, {})
        for aa in translate_dict.values():
            aa_fitness = []
            for codon in aa_to_codons(aa, translate_dict):
                codon_fitness = codon_fitness_scores[pos][codon]
                if codon_fitness is not None:
                    aa_fitness.append(float(codon_fitness[0]))
            if len(aa_fitness) == 0:
                aa_fitness_scores[pos][aa] = None
            else:
                aa_fitness_scores[pos][aa] = float(average(aa_fitness))
    return aa_fitness_scores


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Given either a dictionary of fitness values for barcodes or codons with and without a perturbation, calculates interactions.")
    parser.add_argument('--barcode_fitness', type=str, default=None, nargs=2, metavar=("unperturbed_barcode_fitness", "perturbed_barcode_fitness") help='pickle file encoding a dictionary of barcodes with counts')
    parser.add_argument('--codon_fitness', type=str, default=None, nargs=2, metavar=("unperturbed_codon_fitness", "perturbed_codon_fitness") help='pickle file encoding a dictionary of codons with counts')
    parser.add_argument('--allele_dict', type=str, help='pickle file encoding a dictionary of alleles')
    parser.add_argument('--translate_dict', type=str, help='pickle file encoding a dictionary for translating from codon to amino acid')
    parser.add_argument('--wt_codon_dict', type=str, help='pickle file encoding a dictionary with the wild-type codon for each position')
    parser.add_argument('--wt_time_constants', type=float, nargs=2, metavar=("unperturbed_wt_time_constant, perturbed_wt_time_constant"), help='time constants for wildtype under both conditions')
    parser.add_argument('--aa_index', type=str, help='pickle file encoding a dictionary of amino acids and corresponding indeces used for arranging while plotting')
    args = parser.parse_args()

    

    allele_dict = pickle.load(open(args.allele_dict, 'rb')) 
    translate_dict = pickle.load(open(args.translate_dict, 'rb'))
    wt_codon_dict = pickle.load(open(args.wt_codon_dict, 'rb'))
    unpert_dw, pert_dw = wt_time_constants

    if args.barcode_fitness is not None:
        unpert_barcode_rel_fitness = pickle.load(open(args.barcode_fitness[0], 'rb')) #key is barcode. value is tuple with slope and standard error
        pert_barcode_rel_fitness = pickle.load(open(args.barcode_fitness[1], 'rb')) 
        unpert_codon_rel_fitness = codon_fitness_from_barcodes(unpert_barcode_rel_fitness, allele_dict, wt_codon_dict, translate_dict)
        pert_codon_rel_fitness = codon_fitness_from_barcodes(pert_barcode_rel_fitness, allele_dict, wt_codon_dict, translate_dict)
    else:
        unpert_codon_rel_fitness = pickle.load(open(args.codon_fitness[0], 'rb')) #position->codon->(slope, std_error)
        pert_codon_rel_fitness = pickle.load(open(args.codon_fitness[1], 'rb')) #position->codon->(slope, std_error)

    unpert_aa_rel_fitness = aa_fitness_from_barcodes(unpert_codon_rel_fitness, wt_codon_dict, translate_dict, aa_index)
    pert_aa_rel_fitness = aa_fitness_from_barcodes(pert_codon_rel_fitness, wt_codon_dict, translate_dict, aa_index)


    
    