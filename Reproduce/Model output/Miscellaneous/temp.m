       profile on
        y(:,id) = y(:,id-1) + yp * dt;
        
        % for n, also need to shift
        dx = yp(end-3) * dt;
        xi = parms.xi - dx;
        n = y(1:length(parms.xi),id);
        
        nshift = interp1(parms.xi(:), n, xi(:));
        nshift(isnan(nshift)) = 0;
        nshift(nshift<0) = 0;
        
        y(1:length(parms.xi),id) = nshift;
        
        t2(id) = t2(id-1) + dt;
        
        F2(id) = trapz(parms.xi, parms.xi .* nshift') + trapz(parms.xi, nshift');
        
        profile viewer