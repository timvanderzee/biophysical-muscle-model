function[p, q, Q1] = get_pq(Q0, Q2, Fce)

% mean and standard deviation
% k   = parms.K;
% Q00 = log(1+exp(Q0*k))/k; % note: goes to inf for large k, may need another function
% Q00 = max(Q0, 1e-6);
Q00 = Q0;
Q1 = Fce - Q0;
p = Q1./Q00; 

% limit p
p = 10*tanh(p/10);

q = Q2./Q00 - p.^2;  
% q = log(1+exp(q*k))/k;

q = 1*tanh(q/1); 
q = max(q, 0);

if sum(~isreal(q)) > 0
    keyboard
end
end