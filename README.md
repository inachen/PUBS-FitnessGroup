PUBS-FitnessGroup
=================

##aa_fitness_from_file.py

###Description

Given relative barcode fitnesses, calculates relative codon and amino acid fitnesses, as well as sequence
entropy and information content. Also, given time constants for perturbed and unperturbed wildtype, calculates epistatic interaction of each mutant with the perturbation. Currently propagates variance from standard error of barcode fitnesses, but implements an empty function for estimating variance from a fitness to variance distribution.

Expects barcode fitness to be provided as a dictionary of barcodes matched to a list where the first two values are fitness and standard error, respectively.

Outputs in two different file formats, pickles for downstream visualization and CSVs for hypothesis testing. Non-numbers are replaced by -100 in pickles.

###Example Usage
```bash
python scripts/aa_fitness_from_file.py -h
```

```bash
python scripts/aa_fitness_from_file.py test_files/test_barcode_fitness_1.pkl test_files/test_barcode_fitness_2.pkl --wt_time_constants 0.33 0.2 --allele_dict input_files/allele_dic_with_WT.pkl --translate_dict input_files/translate.pkl --wt_codon_dict input_files/wt_codon_dict.pkl --aa_index input_files/aminotonumber.pkl --weighted_mean --codon_fitness_pickle codon_fitness.pkl --rel_fitness_csv rel_fitness.csv --rel_fitness_pickle rel_fitness.pkl --rel_fitness_variance_csv rel_fitness_variance.csv --sequence_entropy_pickle sequence_entropy.pkl --information_content_pickle information_content.pkl --interaction_pickle interaction.pkl
```

##diff.r

###Description


###Example Usage


##cluster_perturbations.py

###Description

Given amino acid fitness matrices or interaction matrices (pickles) for different perturbations/days as input files, clusters perturbations. First, only matrix items which are real numbers in all matrices are considered, second, distances in "fitness space" or "interaction space" are calculated between each matrix, resulting in a distance matrix. Alternatively, distances can be calculated between mean values for each amino acid position. Finally, perturbations are hierarchically clustered based on distance matrix. The results are visualized as a dendrogram of position mean values.

###Example Usage

```bash
python scripts/cluster_perturbations.py pert1_fitness.pkl pert2_fitness.pkl pert3_fitness.pkl --data_type rel_fitness --pert_names pert1 pert2 pert3 --out_plot clustered_perts_fitness.png
```

```bash
python scripts/cluster_perturbations.py pert1_interaction.pkl pert2_interaction.pkl pert3_interaction.pkl --data_type interaction --pert_names pert1 pert2 pert3 --out_plot clustered_perts_interaction.png
```
