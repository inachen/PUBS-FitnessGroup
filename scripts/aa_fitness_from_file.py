#!/bin/python
import os
import sys
import cPickle as pickle
#from scipy.stats import wald
#import scipy
import numpy as np
import argparse
import pprint

'''
To-do:
    -check on how NaN will be provided
    -allow input and propagation of variance of dw
'''

def array_mean(l, variances=None):
    '''Calculates simple average of list. If variances for each value are provided, also returns combined variance'''
    mean = 0
    var = 0
    n = len(l)
    if n > 0:
        mean = float(sum(l)/n)
        if variances is not None and len(variances) == n:
            var = sum(variances)/n**2
            return mean, var
        else:
            return mean
    return None

def array_weighted_mean(l, variances):
    '''Calculates the inverse variance weighted average of list and returns average with variance'''
    vals = np.array(l)
    real_vals_indices = np.where(~np.isnan(vals))
    real_vals = vals[real_vals_indices]
    variances = np.array(variances)
    real_variances = variances[real_vals_indices]
    n = len(l)
    if n > 0 and n == len(real_variances):
        new_var = float(1/np.sum(1/real_variances**2))
        mean = float(new_var * np.sum(real_vals/real_variances**2))
        return mean, new_var
    return None

def array_sum(mA, mB, vA=None, vB=None):
    """Adds matrices. If matrices with variance for each value provided, returns variance of new matrix.
    (assumes independence of variables)"""
    if vA is None or vB is None:
        return mA+mB
    else:
        return mA+mB,vA+vB

def array_diff(mA, mB, vA=None, vB=None):
    """Subtracts matrices. If matrices with variance for each value provided, returns variance of new matrix.
    (assumes independence of variables)"""
    if vA is None or vB is None:
        return mA-mB
    else:
        return mA-mB,vA+vB

def array_prod(mA, mB, vA=None, vB=None):
    """Multiplies matrices. If matrices with variance for each value provided, returns estimate of variance
    of new matrix. (assumes normality/independence of variables)"""
    if vA is None or vB is None:
        return mA*mB
    else:
        return mA*mB,(mB**2)*vA + (mA**2)*vB

def array_scalar_prod(m, a, v=None):
    """Multiplies scalar by matrix. If matrix with variance for each value provided, returns variance of
    new matrix"""
    if v is None:
        return a*m
    else:
        return a*m,(a**2)*v

def array_ratio(mA, mB, vA=None, vB=None):
    """Divides mA by mB. If matrices with variance for each value provided, returns estimate of variance
    of new matrix. (assumes normality/independence of variables)"""
    if vA is None or vB is None:
        return mA/mB
    else:
        return mA/mB,(mB**2)*vA + (mA**2)*vB

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
    returns as a dictionary of positions, then codons, where the corresponding value is a fitness score'''
    codon_fitnesses = {}

    for pos in wt_codon_dict.keys():
        codon_fitnesses.setdefault(pos,{})
        for codon in translate_dict.keys():
            barcodes = codon_to_barcodes(pos, codon, allele_dict)
            barcode_tuple = [barcode_fitness[x] for x in barcodes if not np.isnan(barcode_fitness[x][0])]
            if len(barcode_tuple) > 0:
                barcode_fitnesses = np.array([x for x,y in barcode_tuple])
                barcode_variances = np.array([y for x,y in barcode_tuple])
                codon_fitness, codon_variance = array_weighted_mean(barcode_fitnesses, barcode_variances)
                codon_fitnesses[pos][codon] = float(codon_fitness)
            else:
                codon_fitnesses[pos][codon] = np.nan
    return codon_fitnesses

def calculate_aa_fitness(codon_fitness, wt_codon_dict, translate_dict, aa_index):
    '''Calculates the fitness of amino acids as the mean fitness of synonymous codons. Returns a numpy array
    of fitness values.'''
    aa_fitnesses = np.empty(shape=(len(wt_codon_dict.keys()),len(aa_index.keys())))

    for pos in wt_codon_dict.keys():
        for aa,i in sorted(aa_index.items(), key=lambda x: x[1]):
            syn_codons = aa_to_codons(aa, translate_dict)
            syn_codon_fitness = np.array([codon_fitness[pos][codon] for codon in syn_codons if codon in codon_fitness[pos]])
            if len(syn_codon_fitness) > 0:
                aa_fitness = np.nanmean(syn_codon_fitness)
            else:
                aa_fitness = np.nan
            aa_fitnesses[pos-1,i] = float(aa_fitness)
    return aa_fitnesses

def wald_test(array, alpha=0.05, temp_file='temp_out.csv'):
    #write out file
    pass
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Given either a dictionary of fitness values for barcodes or codons with and without a perturbation, calculates interactions.")
    parser.add_argument('--barcode_fitness', type=str, default=None, nargs=2, metavar=("unperturbed_barcode_fitness", "perturbed_barcode_fitness"), help='pickle file encoding a dictionary of barcodes with counts')
    parser.add_argument('--codon_fitness', type=str, default=None, nargs=2, metavar=("unperturbed_codon_fitness", "perturbed_codon_fitness"), help='pickle file encoding a dictionary of codons with counts')
    parser.add_argument('--allele_dict', type=str, help='pickle file encoding a dictionary of alleles')
    parser.add_argument('--translate_dict', type=str, help='pickle file encoding a dictionary for translating from codon to amino acid')
    parser.add_argument('--wt_codon_dict', type=str, help='pickle file encoding a dictionary with the wild-type codon for each position')
    parser.add_argument('--wt_time_constants', type=float, nargs=2, metavar=("unperturbed_wt_time_constant, perturbed_wt_time_constant"), help='time constants for wildtype under both conditions')
    parser.add_argument('--aa_index', type=str, help='pickle file encoding a dictionary of amino acids and corresponding indeces used for arranging while plotting')
    args = parser.parse_args()

    allele_dict = pickle.load(open(args.allele_dict, 'rb')) 
    translate_dict = pickle.load(open(args.translate_dict, 'rb'))
    dna_translate_dict = {rna_to_dna(x):y for x,y in translate_dict.items()}
    wt_codon_dict = pickle.load(open(args.wt_codon_dict, 'rb'))
    aa_index = pickle.load(open(args.aa_index, 'rb'))
    unpert_dw, pert_dw = args.wt_time_constants

    if args.barcode_fitness is not None:
        unpert_barcode_rel_fitness = pickle.load(open(args.barcode_fitness[0], 'rb')) #key is barcode. value is tuple with slope and standard error
        pert_barcode_rel_fitness = pickle.load(open(args.barcode_fitness[1], 'rb')) 
        unpert_codon_rel_fitness = codon_fitness_from_barcodes(unpert_barcode_rel_fitness, allele_dict, wt_codon_dict, dna_translate_dict)
        pert_codon_rel_fitness = codon_fitness_from_barcodes(pert_barcode_rel_fitness, allele_dict, wt_codon_dict, dna_translate_dict)
    else:
        unpert_codon_rel_fitness = pickle.load(open(args.codon_fitness[0], 'rb')) #position->codon->(slope, std_error)
        pert_codon_rel_fitness = pickle.load(open(args.codon_fitness[1], 'rb')) #position->codon->(slope, std_error)

    unpert_aa_rel_fitness  = calculate_aa_fitness(unpert_codon_rel_fitness, wt_codon_dict, dna_translate_dict, aa_index)
    pert_aa_rel_fitness = calculate_aa_fitness(pert_codon_rel_fitness, wt_codon_dict, dna_translate_dict, aa_index)

    unpert_aa_norm_fitness = unpert_aa_rel_fitness + 1
    pert_aa_abs_fitness = (pert_aa_rel_fitness+1)*pert_dw
    pert_aa_norm_fitness = pert_aa_abs_fitness/unpert_dw
    #pert_aa_abs_fitness, pert_aa_abs_variance = array_scalar_prod(pert_aa_rel_fitness+1, pert_dw, pert_aa_rel_variance)
    #pert_aa_norm_fitness, pert_aa_norm_variance = array_scalar_prod(pert_aa_abs_fitness, 1/unpert_dw, pert_aa_abs_variance)

    #interaction_product_term, interaction_product_variance = array_prod(pert_aa_norm_fitness, unpert_aa_norm_fitness, pert_aa_norm_variance, unpert_aa_norm_variance)
    #interactions, interaction_variance = array_diff(pert_aa_norm_fitness, interaction_product_term, pert_aa_norm_variance, interaction_product_variance)
    interactions = pert_aa_norm_fitness - pert_aa_norm_fitness*unpert_aa_norm_fitness
    interactions[np.where(np.isnan(interactions))]=-100


