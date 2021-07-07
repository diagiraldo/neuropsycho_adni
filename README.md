# Domain specific composite scores for cognitive impairment

This repository contains:

* Code to calculate the domain composite scores in **new data**, assign a MCI subgroup and predict progression from MCI to dementia within different time windows.
* Code to reproduce the analysis and results presented in the manuscript *"Quantification of cognitive impairment to characterize heterogeneity of patients at risk of developing Alzheimer’s disease dementia"*

### Built With

* R (version 3.6.3)
* Neuropsychological data from the Alzheimer’s Disease Neuroimaging Initiative (ADNI) database. See <http://adni.loni.usc.edu/>

## Reproducing results and plots in the manuscript

### Pre-requisites

To fully reproduce all the analyses you will need the following R libraries:

```R
library(ADNIMERGE)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(lavaan)
library(NbClust)
library(cluster)
library(survival)
library(survminer)
library(randomForest)
library(ROCR)
```

### Folder organization

- Scripts for each of the following steps are in the folder `scripts`.
- Processed ADNI data for each step is saved in a folder named `processed_data`, not included in this public repository but can be generated if the user has access to the `ADNIMERGE` R package.
- Estimated parameters and results are saved in the `results` folder.
- Graphics in EPS and PNG are saved in `plots`.

### Steps

1. Pre-process neuropsychological data in `ADNIMERGE` package with `preprocess_neuropsychological_adni_data.R`
2. Reverse some of the neuropsychological sub-scores with `reverse_scores_neuropsychological_adni_data.R`
3. Pre-process diagnostic data in `ADNIMERGE` package with `preprocess_diagnostic_adni_data.R`
4. Select first visit with complete information and split data for analyses with `select_and_split_neuropsychological_adni_data.R`
5. Calculate Standardized Regression Based (SRB) z-scores with `srb_zscores_neuropsychological_adni_data.R`
6. Perform Confirmatory Factor analysis and calculate domain specific scores with `cfa_with_srb_neuropsychological_adni_data.R`
7. Cluster analysis using domain scores with `mci_clustering_domain_scores.R`, it includes survival analysis to compare progression risk between MCI sub-groups.
8. Automated prediction of progresion from MCI to AD dementia with `mci_progression_predicton_domain_scores.R`.

<!-- LICENSE -->
## License

Distributed under the MIT License.
See `LICENSE` for more information.

<!-- CONTACT -->
## Contact

Diana L. Giraldo - [@diagiraldo](https://github.com/diagiraldo) - <dlgiraldof@unal.edu.co>
