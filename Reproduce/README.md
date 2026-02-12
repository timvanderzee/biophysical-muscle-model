## Reproduce
To reproduce the preprint figures, run `reproduce_manuscript_figures.m`. 

**Table 1: Figure dependencies**
| Figure | GitHub     | Model (time-series) | Data (time-series)    |
|----------|:--------:|:--------:|:--------:|
| Figures 2-4         | ☑️  | ☑️             | Optional |
| Figure 5         | ☑️     | Not required    | Not required |
| Figure 6         | ☑️   | Not required    | Not required |
| Figure 7         | ☑️  | ☑️             | Optional |
| Figure 8         | ☑️  | Not required    | Not required |

Note: Figures 2-4 and Figure 7 require model time-series. These are not provided, but can be obtained from the provided models and parameter values in this repository, explained below. The data time-series is only provided to manuscript reviewers (as supplementary information), and will be available open-access upon publication. Figures 2-4 and Figure 7 can still be produced without the data, but then only the model traces are shown. 

### Obtaining model time-series
To obtain the model time-series (needed for some of the figures, see **Table 1: Figure dependencies**), run `reproduce_model_fits.m`. 

Note: running this script for all fibers, models and methods should result in a folder containing all model time-series. The model time-series should be identical to those from the pre-print, provided that the model code and parameters are not changed. The folder in which time-seris are saved is specified by the user in `outputfolder`. If this folder is then used as `modelfolder` in `reproduce_manuscript_figures.m`, the model traces of Figures 2-4 and Figure 7 can be produced. 

### Evaluate model predictions
If you have obtained the model time-series (see above), you can reproduce the model short-range stiffness (SRS) predictions. If you additionally have the experimental time-series (reviewers only), you can also reproduce the model force RMSDs. Note: SRS model predictions and force RMSDs are already provided in this repository (see `biophysical-muscle-models/Model output`). Reproducing them is just meant to verify that the model time-series calculated above indeed yield the SRS predictions and RMSDs provided in this repository.

### Reproduce short-range stiffness predictions
To reproduce model SRS predictions, run `reproduce_model_SRS.m`

Note: running this script for all fibers, models and methods should replace all the SRS predictions in this repository with new values. These should be identical to the ones provided. 

### Reproduce model force RMSDs
To reproduce model force RMSDs, run `reproduce_model_RMSD.m`

Note: running this script for all fibers, models and methods should replace all the RMSDs in this repository with new values. These should be identical to the ones provided. 
