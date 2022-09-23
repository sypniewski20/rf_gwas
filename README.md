# rf_gwas

1. Clone this repository

2. Start R project in this repo

3. Type `packrat::restore()` for easy setup of required packages


### Usage

1. First run `class_rf_optimise.R` to find best hyperparameters
2. Run `class_run_rf.R` to find best predictor loci (output file `best_predictior_loci.tsv`)

Code based on approach of Brieuc et al. 2018, most changes relates to optimization for multithreading and reproducibility.


### References

Brieuc MSO, Waters CD, Drinan DP, Naish KA. A practical introduction to Random Forest for genetic association studies in ecology and evolution. Mol Ecol Resour. 2018 Jul;18(4):755-766. doi: 10.1111/1755-0998.12773. Epub 2018 Mar 31. PMID: 29504715.
