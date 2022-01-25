t=2*pi/1000:2*pi/1000:20*pi;
x = cos(2*pi*t)+sqrt(-1)*sin(2*pi*t);

dat=[];
dat.label{1}='tmp';
dat.avg = x;
dat.dimord='chan_time';
dat.time=t;

cfg=[];
cfg.taper='hanning';
cfg.output='fourier';
cfg.method='mtmconvol';
cfg.toi = t(1:10:end);
cfg.pad='nextpow2';
cfg.foi=0:0.1:10;
cfg.t_ftimwin=2*pi/10*ones(1,length(cfg.foi));
f=ft_freqanalysis(cfg, dat);

fsel=nearest(f.freq,1);

figure; plot(f.freq, nanmean(squeeze(abs(f.fourierspctrm)),2))
figure; plot(cfg.toi, squeeze(angle(f.fourierspctrm(1,1,fsel,:))))

cfg2=[];
cfg2.taper='hanning';
cfg2.output='pow';
cfg2.method='mtmfft';
f2=ft_freqanalysis(cfg2, f);

figure; plot(f2.freq,f2.powspctrm(fsel,:));

