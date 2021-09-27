# Domain specific composite scores for cognitive impairment

This repository contains:

* Code to calculate the domain composite scores in **new data**, assign a MCI subgroup and predict progression from MCI to dementia within different time windows.
* Code to reproduce the analysis and results presented in the manuscript [*"Quantification of cognitive impairment to characterize heterogeneity of patients at risk of developing Alzheimer’s disease dementia"*](https://alz-journals.onlinelibrary.wiley.com/doi/10.1002/dad2.12237)

### Built With

* R (version 3.6.3)
* Neuropsychological data from the Alzheimer’s Disease Neuroimaging Initiative (ADNI) database. See <http://adni.loni.usc.edu/>

## Calculating scores for new data

### Pre-requisites

* Sub-scores from neuropsychological tests matching the ADNI convention for naming variables. You can input your data to the `example/example_data.csv` file.
* `dplyr` R library should be installed to process data and calculate scores.
* If you want to make predictions about progression from MCI to dementia, you will need the `randomForest` library and the list of pre-trained classifiers in `example/pretrained_RF_MCIprogression_prediction.RData`.

### Run example

In the R console

1. Set the working directory in the repository folder
  ```r
  setwd("~/neuropsycho_adni")
  ```

2. Load the required functions and lists of neuropsychological tests
  ```r
  source("scripts/example_functions.R")
  load("example/tests_lists.RData")

  sel_tests <- c("ADAS", "MMSE", "MOCA", "CLOCK", "TMT", "LM", "CATFL", "AVLT", "BNTMINT")
  it <- unname(unlist(selitems[sel_tests]))
  ```

3. Read the example data and reverse sub-scores that should be reversed
  ```r
  A <- read.csv("example/example_data.csv")
  A <- reverse_scores(A, selitems)
  ```

4. Calculate the Standardized Regression Based (SRB) z-scores using the parameters estimated from a sample of cognitively unimpaired participants in ADNI
  ```r
  srbcoef <- read.csv("results/srb_parameters.csv")
  A_srb <- calculate_srbz(A, it, srbcoef)
  ```

5. Calculate the domain-specific composite scores using the learnt parameters
  ```r
  load("results/cfa_estimates.RData")
  namesfac <- colnames(LW$L)
  S <- composite_scores(newdata = A_srb, weights = LW$W, centerit = meansub)
  S <- cbind(dplyr::select(A_srb, -all_of(it)), S)
  ```

6. Assign one MCI subgroup to each subject in `example_data`. This is done by calculating the distance to each subgroup's representant (medoids) and selecting the closest one.
```r
grmeds <- read.csv("results/medoids_domainscores_MCIsubgroups_k4.csv")
MCIGR <- distance2meds(S, LW$Cz, grmeds, namesfac)
S <- merge(S, select(MCIGR, ID, GR))
```

7. Predict progression from MCI to dementia within the next 1 to 5 years
  ```r
  library(randomForest)
  load("example/pretrained_RF_MCIprogression_prediction.RData")
  PRED <- predict_MCIprogression(S, RFlist, namesfac)
  S <- merge(S, select(PRED, ID, ends_with("prediction")))
  ```

## Reproducing results and plots in the manuscript

### Pre-requisites

To fully reproduce all the analyses you will need the following R libraries:

```r
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
