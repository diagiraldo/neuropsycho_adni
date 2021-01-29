# neuropsycho_adni
Code for the analysis of neuropsychological data in ADNI

## List of folders
- scripts
- ADNI original data (not here)
- processed data and results (not here)
- plots

## Steps

1. Pre-process neuropsychological data in `ADNIMERGE` package with `preprocess_neuropsychological_adni_data.R`
2. Reverse some of the neuropsychological sub-scores with `reverse_scores_neuropsychological_adni_data.R`
3. Pre-process diagnostic data in `ADNIMERGE` package with `preprocess_diagnostic_adni_data.R`
4. Select first visit with complete information and split data for analyses with `select_and_split_neuropsychological_adni_data.R`
5. Calculate Standardized Regression Based (SRB) z-scores with `srb_zscores_neuropsychological_adni_data.R`
6. Perform Confirmatory Factor analysis and calculate domain specific scores with `cfa_with_srb_neuropsychological_adni_data.R`