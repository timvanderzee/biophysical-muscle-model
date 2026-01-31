# biophysical-muscle-model repository

Welcome to this repository! 

This repository allows you to use a range of biophysical muscle models and evaluate their predictions of short-range stiffness against experimental data. The models and data are described in detail in our [preprint](https://www.biorxiv.org/content/10.1101/2025.10.31.685881v1). 

There are three types of things you could use this code for:
1. Reproduce the [preprint's](https://www.biorxiv.org/content/10.1101/2025.10.31.685881v1) figures.
2. Evaluate model predictions
3. Refit model parameters
4. Adapt the model(s) to fit your own needs

We recommend that you follow the order shown above. Next, we will provide detail on each of these steps.

## 0. Set up paths
Before you get started, you need to set up some paths. Please edit the following lines of code in `get_paths.m`:

```bash
% this datafolder should contain the subfolders 2017 and 2018
datafolder = '';

% this folder is where the model fits are or will be saved
modelfolder = '';

% this folder should contain the biophysical-muscle-model repository
githubfolder = '';
```

Note: not all figures require all three folders (see Table below)

**Table 1: Figure dependencies**
| Figure | GitHub     | Model (time-series) | Data (time-series)    |
|----------|:--------:|:--------:|:--------:|
| Figures 2-4         | ☑️  | ☑️             | Optional |
| Figure 5         | ☑️     | Not required    | Not required |
| Figure 6         | ☑️   | Not required    | Not required |
| Figure 7         | ☑️  | ☑️             | Optional |
| Figure 8         | ☑️  | Not required    | Not required |

## 1. Reproduce the preprint figures
To reproduce the preprint figures, run `reproduce_manuscript_figures.m`. 

Note: Figures 2-4 and Figure 7 require model time-series. These are not provided, but can be obtained from the provided models and parameter values in this repository, explained below. The data time-series is only provided to manuscript reviewers (as supplementary information), and will be available open-access upon publication. Figures 2-4 and Figure 7 can still be produced without the data, but then only the model traces are shown. 

### Obtain model predictions of time-series
To obtain the model time-series, run `reproduce_model_fits.m`. 

Note: running this script for all fibers, models and methods should result in a folder containing all model time-series. The model time-series should be identical to those from the pre-print, provided that the model dynamics and model parameters are not changed. The folder in which time-seris are saved is specified by the user in `outputfolder`. If this folder is then used as `modelfolder` in `reproduce_manuscript_figures.m`, the model traces of Figures 2-4 and Figure 7 can be produced. 

## 2. Evaluate model predictions
If you have evaluated the models and obtained the model predictions of time-series (see above), you can reproduce the short-range stiffness (SRS) model predictions. If you additionally have the experimental time-series (reviewers only), you can also reproduce the model force RMSDs. Note: SRS model predictions and force RMSDs are already provided in this repository (i.e. in `biophysical-muscle-models/Model output`). Reproducing them is just meant to verify that the model time-series calculated above indeed yield the SRS precitions and RMSDs provided in this repository.

### Reproduce short-range stiffness predictions
To reproduce model SRS predictions, run `reproduce_model_SRS.m`

Note: running this script for all fibers, models and methods should replace all the SRS predictions in this repository with new values. These should be identical to the ones provided. 

### Reproduce model force RMSDs
To reproduce model force RMSDs, run `reproduce_model_RMSD.m`

Note: running this script for all fibers, models and methods should replace all the RMSDs in this repository with new values. These should be identical to the ones provided. 

## 3. Refit model parameters


## 4. Adapt the model(s) to fit your own needs
