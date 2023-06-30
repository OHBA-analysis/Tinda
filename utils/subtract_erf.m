function data = subtract_erf(data, erf, events, twindow)

if ~exist('twindow', 'var')  || isempty(twindow)
    twindow = [0 size(erf,2)-1];
end

if length(size(erf)) == 3
    erf = nanmean(erf,3);
end

for e=1:size(events,1)
    data(:, events(e,1)+twindow(1) : events(e,1)+twindow(2)) = data(:, events(e,1)+twindow(1) : events(e,1)+twindow(2)) - erf;
end