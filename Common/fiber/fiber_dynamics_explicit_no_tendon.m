function[dx] = fiber_dynamics_explicit_no_tendon(t,y, parms)

% Get velocity and calcium
if numel(parms.vts) == 1
    vMtilda = parms.vts;
else
    vMtilda = interp1(parms.ti, parms.vts, t);
end

if numel(parms.Cas) == 1
    Ca = parms.Cas;
else
    Ca = interp1(parms.ti, parms.Cas, t);
end

% activation
Act = parms.actfunc(Ca, parms);

% States
Q0  = y(1);
Q1  = y(2);
Q2  = y(3);
L  = y(4);
Non = y(5);
DRX = y(6);
R = y(7);

% Cross-bridge states
k   = parms.K;
F   = Q1 + Q0;
F   = log(1+exp(F*k))/k;
Q00 = log(1+exp(Q0*k))/k;

% mean and standard deviation
p = Q1./Q00; 
q = Q2./Q00 - p.^2;  
q = log(1+exp(q*k))/k;

% Thin and thick filament
if (parms.kon == 0) && (parms.koff == 0) && (parms.koop == 0)
    dNondt = ((Act - Non) / .005);
else 
    [Jon, Joff] = ThinFilament_Dynamics(Act, Q0, Non, parms.kon, parms.koff, parms.koop, 1);
    dNondt = Jon - Joff;
end

[J1, J2] = ThickFilament_Dynamics(Q0, F, DRX, parms.J1, parms.J2, parms.JF, 1, R);

% without limit
% IGef{1} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2)));
% IGef{2} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2))).*(c(2,:)-c(3,:)*k(2)/2);
% IGef{3} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2))).*((c(2,:)-c(3,:)*k(2)/2).^2+c(3,:)/2);

% with limit
IGef{1} = @(c,k)(c(1,:)*k(1).*exp(10*tanh((c(3,:)*k(2)^2/4-c(2,:)*k(2))/10)));
IGef{2} = @(c,k)(c(1,:)*k(1).*exp(10*tanh((c(3,:)*k(2)^2/4-c(2,:)*k(2))/10))).*(c(2,:)-c(3,:)*k(2)/2);
IGef{3} = @(c,k)(c(1,:)*k(1).*exp(10*tanh((c(3,:)*k(2)^2/4-c(2,:)*k(2))/10))).*((c(2,:)-c(3,:)*k(2)/2).^2+c(3,:)/2);

if parms.approx
    % erf approximation
    erfap = @(x) (exp(x)-exp(-x)) ./ (exp(x)+exp(-x));
    IG{1} =  @(x,c)1/2*c(1,:).*erfap((x-c(2,:))./sqrt(c(3,:)));
    IG{2} =  @(x,c)1/2*c(1,:).*erfap((x-c(2,:))./sqrt(c(3,:))) .* c(2,:)-1/2.*sqrt(c(3,:)).*c(1,:)/sqrt(pi).*exp(-(x-c(2,:)).^2./c(3,:));
    IG{3} =  @(x,c)1/2*c(1,:).*erfap((x-c(2,:))./sqrt(c(3,:))) .* (c(2,:).^2+c(3,:)/2)-1/2*sqrt(c(3,:)).*c(1,:)/sqrt(pi).*exp(-(x-c(2,:)).^2./c(3,:)).*(c(2,:)+min(x,1e4));

else
    % erf function
    IG{1} =  @(x,c)1/2*c(1,:).*erf((x-c(2,:))./sqrt(c(3,:)));
    IG{2} =  @(x,c)1/2*c(1,:).*erf((x-c(2,:))./sqrt(c(3,:))).*c(2,:)-1/2.*sqrt(c(3,:)).*c(1,:)/sqrt(pi).*exp(-(x-c(2,:)).^2./c(3,:));
    IG{3} =  @(x,c)1/2*c(1,:).*erf((x-c(2,:))./sqrt(c(3,:))).*(c(2,:).^2+c(3,:)/2)-1/2*sqrt(c(3,:)).*c(1,:)/sqrt(pi).*exp(-(x-c(2,:)).^2./c(3,:)).*(c(2,:)+min(x,1e4));

end

% points where integrals is evaluated
k1 = [parms.k11 parms.k12];
k2 = [parms.k21 -parms.k22];

% Compute Qdot
[Q0dot, Q1dot, Q2dot, Rdot] = CrossBridge_Dynamics(Q0, p, q, parms.f, parms.w, k1, k2, IGef, Non, DRX, IG, parms.b, parms.k, R, parms.dLcrit, parms.ps2);

% velocity - independent derivative
F0dot  = Q1dot + Q0dot;

kse = parms.kse * (F + parms.kse0);
Ld  = (vMtilda .* parms.gamma .* kse - F0dot) ./ (Q0 + kse);

dQ0dt = Q0dot;
dQ1dt = (Q1dot + 1 * Ld .* Q0);
dQ2dt = (Q2dot + 2 * Ld .* Q1);

% note: the term "Q0dot + Rdot" is the part of Q0dot due to regular
% attachment. explanation: if Rdot is great, it means that Q0 is losing to
% R. this implies that for a given Q0dot, the part due to regular
% attachment must be greater. thus, we need to add Rdot to Q0dot
dDRXdt = J1 - J2 - (Q0dot + Rdot);

dx = [dQ0dt; dQ1dt; dQ2dt; Ld; dNondt; dDRXdt; Rdot];

end