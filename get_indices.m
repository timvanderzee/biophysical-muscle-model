function[id0,id1,id2] = get_indices(t, tiso, ts, dTt, dTc, ISI, Ca)

dt = median(diff(t));
N = round(.01/dt); % 10 ms

% SRS indices
id0 = nan(length(Ca), 10*N);
id1 = nan(length(Ca), N);
id2 = nan(length(Ca), N);

for i = 1:length(Ca)
    if sum(t < (ts(i) + tiso - 3*dTt - 2*dTc - ISI)) > 0
        id0(i,:) = find(t < (ts(i) + tiso - 3*dTt - 2*dTc - ISI),10*N, 'last');
    end
    
    if nargout > 1
        id1(i,:) = find(t > (ts(i) + tiso - 3*dTt),N, 'first');
        id2(i,:) = find(t > (ts(i) + tiso - 3*dTt - 2*dTc - ISI),N, 'first');
    end
end
end