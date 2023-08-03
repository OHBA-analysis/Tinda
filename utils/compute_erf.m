function erf = compute_erf(data, events, twindow, baseline)

%% input 
% data (chan x time)
% events (n_events x [sample event_id])
% twindow [first_sample, last_sample] compared to trigger sample
% baseline [first_sample, last_sample] compared to trigger sample


erf = zeros(size(data,1), diff(twindow)+1, length(events)); % chan x time x events

for e=1:size(events,1)
    erf(:,:,e) = data(:, events(e,1)+twindow(1) : events(e,1)+twindow(2));
end


if exist('baseline', 'var') && ~isempty(baseline)
    bl = zeros(size(data,1), diff(baseline)+1, length(events)); % chan x time x events
    for e=1:size(events,1)
        bl(:,:,e) = data(:, events(e,1)+baseline(1) : events(e,1)+baseline(2));
    end
    erf = erf - nanmean(nanmean(bl,3),2);
end

