function[modelfunc, modelname] = look_up_model(model)

if contains(model, 'XB')
    modelfunc = @fiber_dynamics_explicit_length_v2;
else
    modelfunc = @hill;
end

if strcmp(model, '3-state XB coop')
    modelname = 'biophysical_full_regular';

elseif strcmp(model, '2-state XB coop')
    modelname = 'biophysical_thin_regular';

elseif strcmp(model, '2-state XB')
    modelname = 'biophysical_no_regular';

elseif strcmp(model, '4-state XB coop')
    modelname = 'biophysical_full_alternative';
end
end