#!/bin/python
import os
import sys
import cPickle as pickle
from scipy.stats import wald
import scipy
import numpy as np
import argparse


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
    variances = np.array(variances)
    n = len(l)
    if n > 0 and n == len(variances):
        new_var = float(1/np.sum(1/variances**2))
        mean = float(new_var * np.sum(vals/variances**2))
        return new_var*mean
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
    returns as a dictionary of positions, then codons, where the corresponding value is a tuple with fitness
    score and variance.'''
    codon_fitness_scores = {}

    for pos in wt_codon_dict.keys():
        codon_fitness_scores.setdefault(pos,{})
        for codon in translate_dict.keys():
            barcode_fitnesses = np.array(sorted([x for x in codon_to_barcodes(pos, codon, allele_dict) if not np.isnan(x)]))
            codon_fitness[pos][codon] = float(np.nanmean(barcode_fitnesses))

    return codon_fitness_scores

def calculate_aa_fitness(codon_fitness, wt_codon_dict, translate_dict, aa_index):
    '''Calculates the fitness of amino acids as the mean fitness of synonymous codons. Returns a numpy array
    of fitness values and an array of variances.'''
    aa_fitnesses = np.empty(size=(len(wt_codon_dict),len(translate_dict)))
    aa_variances = np.empty(size=(len(wt_codon_dict),len(translate_dict)))

    for pos in wt_codon_dict.keys():
        for aa,i in sorted(aa_index.items(), key=lambda x: x[1]):
            syn_codons = aa_to_codons(aa, translate_dict)
            syn_codon_fitness = np.array([codon_fitness_scores[pos][codon][0] for codon in syn_codons if codon in codon_fitness_scores[pos]])
            syn_codon_variance = np.array([codon_fitness_scores[pos][codon][1] for codon in syn_codons if codon in codon_fitness_scores[pos]])
            if len(syn_codon_fitness) > 0:
                aa_fitness, aa_variance = array_mean(syn_codon_fitness, syn_codon_variance)
            else:
                aa_fitness, aa_variance = np.nan, np.nan
            aa_fitnesses[pos-1,i] = float(aa_fitness)
            aa_variances[pos-1,i] = float(aa_variance)
    return aa_fitnesses


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
    aa_index = pickle.load(open(args.aa_index, 'rb'))
    unpert_dw, pert_dw = wt_time_constants

    if args.barcode_fitness is not None:
        unpert_barcode_rel_fitness = pickle.load(open(args.barcode_fitness[0], 'rb')) #key is barcode. value is tuple with slope and standard error
        pert_barcode_rel_fitness = pickle.load(open(args.barcode_fitness[1], 'rb')) 
        unpert_codon_rel_fitness = codon_fitness_from_barcodes(unpert_barcode_rel_fitness, allele_dict, wt_codon_dict, translate_dict)
        pert_codon_rel_fitness = codon_fitness_from_barcodes(pert_barcode_rel_fitness, allele_dict, wt_codon_dict, translate_dict)
    else:
        unpert_codon_rel_fitness = pickle.load(open(args.codon_fitness[0], 'rb')) #position->codon->(slope, std_error)
        pert_codon_rel_fitness = pickle.load(open(args.codon_fitness[1], 'rb')) #position->codon->(slope, std_error)

    unpert_aa_rel_fitness,unpert_aa_rel_variance  = aa_fitness_from_barcodes(unpert_codon_rel_fitness, wt_codon_dict, translate_dict, aa_index)
    pert_aa_rel_fitness,pert_aa_rel_variance = aa_fitness_from_barcodes(pert_codon_rel_fitness, wt_codon_dict, translate_dict, aa_index)

    unpert_aa_norm_fitness = unpert_aa_rel_fitness + 1
    unpert_aa_norm_variance = unpert_aa_rel_variance
    pert_aa_abs_fitness, pert_aa_abs_variance = array_scalar_prod(pert_aa_rel_fitness+1, pert_dw, pert_aa_rel_variance)
    pert_aa_norm_fitness, pert_aa_norm_variance = array_scalar_prod(pert_aa_abs_fitness, 1/unpert_dw, pert_aa_abs_variance)

    interaction_product_term, interaction_product_variance = array_prod(pert_aa_norm_fitness, unpert_aa_norm_fitness, pert_aa_norm_variance, unpert_aa_norm_variance)
    interactions, interaction_variance = array_diff(pert_aa_norm_fitness, interaction_product_term, pert_aa_norm_variance, interaction_product_variance)

