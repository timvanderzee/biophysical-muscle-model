# biophysical-muscle-model repository

Welcome to this repository! 

This repository allows you to use a range of biophysical muscle models and evaluate their predictions of short-range stiffness against experimental data. The models and data are described in detail in our [preprint](https://www.biorxiv.org/content/10.1101/2025.10.31.685881v1). 

There are three types of things you could use this code for:
1. Reproduce the [preprint's](https://www.biorxiv.org/content/10.1101/2025.10.31.685881v1) figures.
2. Evaluate model predictions for a range of conditions
3. Adapt the model(s) to fit your own needs

We recommend that you follow the order shown above. Next, we will provide detail on each of these steps.

## 1. Reproduce the preprint figures

To reproduce the preprint figures, run `reproduce_manuscript_figures.m`. 

For this script to succesfully reproduce all figures, it needs:
1. Model predictions. There are two types of model predictions:
    - Time-series (e.g., forces, lengths). These are **not** provided in this repository, but they can be reproduced from models and model parameters (which are provided). Instructions on how to do that are provided later in this README. 
    - Short-range stiffness predictions. There are provided in this repository.
    - Model computational time estimates.
2. The experimental data. There are two types of experimental data:
   - Time-series (e.g., forces, lengths). These are currently **not** available. They will become available upon publication. In the mean time, only manuscript reviewers will have this data, attached as supplementary information to the submitted manuscript. 
   - Short-range stiffness data. This data is included in this repository.
3. Estimates of deviations between model predictions and experimental data. The are provided in this repository.

Here is a quick overview of what each of the preprint's result figures needs:
| Figure | Only requires this repo   | Other requirements     |
|----------|:--------:|:--------:|
| Figures 2-4         | ❌ | Model and data time-series|
| Figure 5         | ☑️ | |
| Figure 6         | ☑️ | |
| Figure 7         | ❌ | Model and data time-series|
| Figure 8         | ☑️| |

### Obtain model predictions of time-series
If this is your first time using this repository, running the `reproduce_manuscript_figures.m` script resulted in reproducing Figures 5, 6 and 8, but not Figures 2-4 or Figure 7 (see Table shown above). The latter figures cannot be fully reproduced at this time, because the experimental time-series are only provided upon publication. However, you can reproduce the model predictions shown in these figures. To accomplish this, you need to evaluate the models. Instructions are provided below:







