function[p, q, Q1] = get_pq(Q0, Q2, Fce)

% compute first order moment
Q1 = Fce - Q0;

% mean and standard deviation of the XB distribution
p = Q1./Q0;                % mean
q = Q2./Q0 - p.^2;         % standard deviation

% constraining p and q
p = 10*tanh(p/10);          % constraining minimum and maximum
q = log(1+exp(q*10))/10;    % constraining to be larger than zero
q = 2*tanh(q/2);            % constraining the maximum 

end