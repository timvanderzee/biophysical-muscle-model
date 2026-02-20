# Reproduce
You can reproduce the following things:
- model parameter fitting
- model force time-series
- model force root-mean-square deviations (RMSDs) from data
- model short-range stiffness (SRS)
- manuscript figures

Each of these are explained below

## Model parameter fitting
To perform parameter fitting you need CasADi software, see: https://web.casadi.org/get/. For more information see [biophysical-muscle-model/Fitting](https://github.com/timvanderzee/biophysical-muscle-model/tree/main/Fitting). 

To reproduce parameter fitting, run the `reproduce_model_fitting.m` script. Running this script for all fibers and models should result in the same parameter values as in our pre-print. If you're just interested in the parameter values, they are provided in the [Parameters folder](https://github.com/timvanderzee/biophysical-muscle-model/tree/main/Reproduce/Parameters).

## Model force time-series
To reproduce model forces, run `reproduce_model_forces.m`. Running this script for all fibers, models and methods should result in a folder containing all model time-series. The model force time-series should be identical to those from the pre-print, provided that the model code and parameters are not changed. The folder in which force time-seris are saved is specified by the user in `outputfolder`.

## Short-range stiffness predictions
To reproduce model SRS predictions, run `reproduce_model_SRS.m`. If you're just interested in the model SRS values, they are provided in the [Model output folder](https://github.com/timvanderzee/biophysical-muscle-model/tree/main/Reproduce/Model%20output).

## Model force root-mean-square deviations
To reproduce model force RMSDs, run `reproduce_model_RMSD.m`. If you're just interested in the force RMSDs values, they are provided in the [Model output folder](https://github.com/timvanderzee/biophysical-muscle-model/tree/main/Reproduce/Model%20output).

## Manuscript figures
To reproduce the preprint figures, run `reproduce_manuscript_figures.m`. Figures 2-4 and Figure 7 require model and data time-series (see Table 1). Model time-series can be obtained from the provided models and parameter values in this repository (see "Reproduce model forces"). The data time-series is only provided to manuscript reviewers (as supplementary information), and will be available open-access upon publication. Figures 2-4 and Figure 7 can still be produced without the data, but then only the model traces are shown. The other figures require model SRS and RMSDs, which are provided in this repository in [Model output folder](https://github.com/timvanderzee/biophysical-muscle-model/tree/main/Reproduce/Model%20output).

**Table 1: Figure dependencies**
| Figure | GitHub     | Model (time-series) | Data (time-series)    |
|----------|:--------:|:--------:|:--------:|
| Figures 2-4         | ☑️  | ☑️             | Optional |
| Figure 5         | ☑️     | Not required    | Not required |
| Figure 6         | ☑️   | Not required    | Not required |
| Figure 7         | ☑️  | ☑️             | Optional |
| Figure 8         | ☑️  | Not required    | Not required |



