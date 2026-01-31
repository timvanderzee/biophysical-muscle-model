# biophysical-muscle-model repository

Welcome to this repository! 

This repository allows you to use a range of biophysical muscle models and evaluate their predictions of short-range stiffness against experimental data. The models and data are described in detail in our [preprint](https://www.biorxiv.org/content/10.1101/2025.10.31.685881v1). 

There are three types of things you could use this code for:
1. Reproduce the [preprint's](https://www.biorxiv.org/content/10.1101/2025.10.31.685881v1) figures.
2. Evaluate model predictions for a range of conditions
3. Refit model parameters
4. Adapt the model(s) to fit your own needs

We recommend that you follow the order shown above. Next, we will provide detail on each of these steps.

## 1. Reproduce the preprint figures

To reproduce the preprint figures, run `reproduce_manuscript_figures.m`. 

In order to run this script without errors, you need to modify the following code: 

```bash
% this datafolder should contain the subfolders 2017 and 2018
datafolder = '';

% this folder should contain model fits
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

Thus, Figures 2-4 and Figure 7 require model time-series. These are not provided, but can be obtained from the provided models and parameter values in this repository, explained below. The data time-series is only provided to manuscript reviewers (as supplementary information), and will be available open-access upon publication. Figures 2-4 and Figure 7 can still be produced without the data, but then only the model traces are shown. 

### Obtain model predictions of time-series
To obtain the model time-series, run `reproduce_model_fits.m`. 

In order to run this script without errors, you need to modify the following code: 

```bash
% in this folder the forces will be saved (user's choice)
outputfolder = '';

% this folder should contain the biophysical-muscle-model repository
githubfolder = '';
```

Running this script for all fibers, models and methods should result in a folder containing all model time-series. The folder is specified by the user in `outputfolder`. If this folder is then used as `modelfolder` in `reproduce_manuscript_figures.m`, Figures 2-4 and Figure 7 can be reproduced. 

## 2. Evaluate model predictions for a range of conditions


## 3. Refit model parameters


## 4. Adapt the model(s) to fit your own needs
