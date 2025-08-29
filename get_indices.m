function[id0,id1,id2] = get_indices(t, ts, dTt, dTc, ISI, Ca)

dt = mean(diff(t));
N = round(.01/dt); % 10 ms

% SRS indices
id0 = nan(length(Ca), 10*N);
id1 = nan(length(Ca), N);
id2 = nan(length(Ca), N);

for i = 1:length(Ca)
    id0(i,:) = find(t < (ts(i) + mean(diff(ts)) - 3*dTt - 2*dTc - ISI),10*N, 'last');
    id1(i,:) = find(t > (ts(i) + mean(diff(ts)) - 3*dTt),N, 'first');
    id2(i,:) = find(t > (ts(i) + mean(diff(ts)) - 3*dTt - 2*dTc - ISI),N, 'first');
end
end