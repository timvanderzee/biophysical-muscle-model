# Fitting

Here, we will show you how to fit parameters on a given force trajectory. 

For fitting, you will need CasADi software, see: https://web.casadi.org/get/. Make sure that the folder containing the downloaded CasADi software is added to your MATLAB path. 

Once you have downloaded CasADi, you can run `fit_any_force_trajectory.m`. This script fits a few parameters on an example force trajectory in 7 steps:
1. specify which model to fit
2. specify which parameters to fit
3. obtain initial parameter values
4. specify which data to fit
5. obtain initial guess for model states
6. do fitting
7. validate fitting through running a simulation with fitted parameters

These steps are explained in further detail below.

## Step 1: specify which model to fit
You can choose between the following models:
- 2-state XB
- 2-state XB coop
- 3-state XB coop
- 4-state XB coop

## Step 2: specify which parameters to fit
As an example, four parameters are selected:
- `f`: the crossbridge attachment rate
- `k11`: the crossbridge detachment rate at negative strain
- `J1`: the super-relaxed state forward rate constant
- `kse`: an SE stiffness parameter

You can fit any of the parameters stored in the `parms` struct. For each of the selected parameters, you need to specify lower and upper bounds in the `bnds` struct. 

## Step 3: obtain initial parameter values
You need to obtain parameter values, including initial values for the parameters selected for fitting. Here, we load in parameters from our pre-print for a selected fiber. 

## Step 4: specify which data we want to fit
Similar to as in `test_model.m`, we need an input struct called `input` that contains the following fields:
-     t: [1×N double]
-     L: [1×N double]
-     v: [1×N double]
-     Ca: [1×N double]

In addition to these fields, we need the following field:
-     F: [1×N double]

This last field contains the force that we want to fit. Here, as an example, we will fit an adjusted simulated force trajectory. We start with the simulated force corresponding to the above input (`input`), the model (`model`), and the model parameters (`oldparms`). Next, we adjust this force through multiplying it with a scaling factor, and adding an offset and noise. As a result of these adjustments, there will be a mismatch between the initial model force and the target model force. This mismatch will then be reduced by re-fitting the model parameters.

## Step 5: obtain initial guess for model states
We need an initial guess for the model states, given the selected model, parameters and input. This initial guess is stored in the variable `IG`. 

## Step 6: do fitting
Here, we can the function `fit_model_parameters_v2`, with the following inputs:
- `model`: model to be fitted
- `oldparms`: parameter values
- `optparms`: parameters that will be fitted
- `bnds`: bounds of the parameters that will be fitted
- `input`: input to the model (calcium, length, velocity, target force)
- `IG`: initial guess of model states
- `weights`: weights for terms in the cost function

 This functions returns the following:
 - `newparms`: the new parameters, including the ones fitted
 - `out`: simulation output (e.g. model states)

## Step 7: validate through simulating with obtained parameters
This is to check that everything went well. 
