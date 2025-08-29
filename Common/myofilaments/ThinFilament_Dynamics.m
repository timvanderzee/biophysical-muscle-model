function[Jon, Joff] = ThinFilament_Dynamics(Ca, Q0, Non, kon, koff, koop, Ntot)

    Noff = Ntot - Non;
    Jon     = Ca .* kon .* Noff  .* (1 + koop .* (Non/Ntot)); % Eq (1)
    Joff    = koff * (Non - Q0) .*     (1 + koop * Noff/Ntot); % Eq (2)
    
end