function [spctrm, freq, time] = fourieranalysis_circleVpath(vpathcircle, trialonset)

% put the unit circle equivalent of the vpath in a FT data struture
fsample=250;
dat = [];
dat.time = (1:length(vpathcircle))/fsample;
dat.dimord = 'chan_time';
for k=1:size(vpathcircle,2);
    dat.label{k,1} = sprintf('vpcircle_long_%d',k);
end
dat.avg = vpathcircle';

% 20 second sliding window, sliding with 10 ms
timres=0.01;
fres=0.05;
cfg=[];
cfg.toi = timres:timres:length(vpathcircle)/fsample;
cfg.taper = 'hanning';
cfg.output = 'fourier';
cfg.foi = fres:fres:5;
cfg.t_ftimwin = (1./fres)*ones(size(cfg.foi));
cfg.method = 'mtmconvol';
f=ft_freqanalysis(cfg, dat);
freq=f.freq;

% downsample trial onset time course and make an evoked spectrum
if exist('trialonset', 'var') && ~isempty(trialonset)
    onsetidx = find(trialonset==1);
else
    onsetidx = fsample+1:fsample:length(f.time)-fsample;
end
W=1;
time = -W:timres:W;
for k=1:length(onsetidx)
    idx(k) = nearest(f.time, dat.time(onsetidx(k)));
end
idx(idx<=W/timres)=[];
spctrm = nan(length(idx), length(dat.label), length(f.freq), length(time));
for k=1:length(idx)
        spctrm(k,:,:,:) = squeeze(f.fourierspctrm(1,:,:,idx(k)-W/timres:idx(k)+W/timres));
end
