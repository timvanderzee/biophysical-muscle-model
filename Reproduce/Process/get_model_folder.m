function[fullmodelfolder, subfolder, modelname] = get_model_folder(modelfolder, mcode, discretized_model)

% folders
% save folder
if mcode(1) == 1 % biophysical model
    if discretized_model 
        subfolder =  'Discretized';
    else
        subfolder = 'Approximated';
    end
else % Hill-type model
    subfolder = 'Hill';
end

% models
models = {'biophysical','Hill', 'PE'};

% model variations
mvars = {'regular','alternative'};

% type of cooperativity
cvars = {'full','thin','no'};

% choose the model version
model = models{mcode(1)};
mvar = mvars{mcode(2)};
cvar = cvars{mcode(3)};

if strcmp(model,'biophysical')
    modelname = [model,'_',cvar,'_',mvar];    
else
    modelname = [model,'_',mvar];
end

% folder containing the .mat files
fullmodelfolder = fullfile(modelfolder, subfolder);

end