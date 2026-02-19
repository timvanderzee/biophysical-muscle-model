# Fitting

Here, we will show you how to fit parameters on a given force trajectory. 
For fitting, you will need CasADi software, see: https://web.casadi.org/get/
Make sure that the folder containing the downloaded CasADi software is added to your MATLAB path. 

Once you have downloaded CasADi, you can run the `fit_any_force_trajectory.m` script.

This script fits a few parameters on an example force trajectory in 7 steps:
1. specify which model to fit
2. specify which parameters to fit
3. obtain initial parameter values
4. specify which data to fit
5. obtain initial guess for model states
6. do fitting
7. validate fitting through running a simulation with fitted parameters

These steps are explained in further detail below

## Step 1: specify which model to fit
You can choose any of the biophysical models in this repository (see biophysical-muscle-model/Test for more info)

## Step 2: specify which parameters to fit
As an example, four parameters are selected:
- f: the crossbridge attachment rate
- k11: the crossbridge detachment rate at negative strain
- J1: the super-relaxed state forward rate constant
- kse: an SE stiffness parameter

You can fit any of the parameters, stored in the struct
