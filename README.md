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

This script takes in two fitness matrices along with the variance matrix and outputs the difference matrix, difference array along position, Wald test p value matrix, and Wald test p value array along position.

Wald test: For each allele (pos, aa), we test the difference of the fitness values for this allele in the two matrices (f<sub>i</sub> - f<sub>j</sub>). The null hypothesis is f<sub>i</sub> - f<sub>j</sub> = 0. A chi-squared distribution is used to compute the p value.

### Files
**Make sure the following directories are present:**  
data_files/ : input csv files containing matrices of fitness values  
out_files/ : where the output csv files will be stored  

(These directories can be changed at the top of the diff.r script.)

###Constants:
(Also located at the top of the diff.r script  
  
File separation, change based on operating system  

```
FSEP = "/"
```

Change significance level (alpha) as needed  
```
ALPHA = 0.05
```

###Usage
####Input
If fitness matrices to be compared is stored in `data_files/fitness1.csv` and `data_files/fitness2.csv`, and the variance matrix is stroed in `data_files/variance.csv`, then

```
run("fitness1", "fitness2", "variance")
```


####Output:  
```
out_files/diff_mat
out_files/diff_lst
out_files/wald_mat
out_files/wald_lst
```



##cluster_perturbations.py

###Description

Given amino acid fitness matrices (pickles) for different perturbations/days as input files, clusters perturbations. First, only matrix items which are real numbers in all matrices are considered, second, distances in "fitness space" are calculated between each matrix, resulting in a distance matrix. Alternatively, distances can be calculated between mean fitnesses for each amino acid position. Finally, perturbations are hierarchically clustered based on distance matrix. The results are visualized as a dendrogram of mean position fitnesses.

###Example Usage
```bash
python scripts/cluster_perturbations.py pert1_fitness.pkl pert2_fitness.pkl pert3_fitness.pkl --pert_names pert1 pert2 pert3
```
