# biophysical-muscle-model repository

Welcome to this repository! 

This repository allows you to use a range of biophysical muscle models and evaluate their predictions of short-range stiffness against experimental data. The models and data are described in detail in our [preprint](https://www.biorxiv.org/content/10.1101/2025.10.31.685881v1). 

There are a few things you could use this code for:
1. To reproduce the [preprint's](https://www.biorxiv.org/content/10.1101/2025.10.31.685881v1) results and figures, go to [Reproduce results](https://github.com/timvanderzee/biophysical-muscle-model/reproduce). 
2. To refit model parameters, go to: 
3. To adapt the model(s) to fit your own needs, see: 

Next, we will provide detail on each of these steps.

## Before you start: set up your local paths
Before you get started, you need to set up some paths. Please edit the following lines of code in `get_paths.m`:

```bash
% this datafolder should contain the subfolders 2017 and 2018
datafolder = '';

% this folder is where the model fits are or will be saved
modelfolder = '';

% this folder should contain the biophysical-muscle-model repository
githubfolder = '';
```

Note: `githubfolder` is always required, `modelfolder` is only required if you are looking to reproduce model time-series, and `datafolder` is only required if you are looking to reproduce data time-series. 

## Contact
If you have any questions regarding this repository, please contact me at tim.vanderzee@kuleuven.be
