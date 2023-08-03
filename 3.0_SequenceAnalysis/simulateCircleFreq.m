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




%%



dat=[];
dat.label{1} = 'tmp';
dat.dimord = '{rpt}_chan_time';
dat.fsample=250;
for k=1:55
  vpos=nan(length(vpath{k}),1);
for ik=1:12
  vpos(vpath{k}==ik) = disttoplot_manual(ik);
end
  dat.trial{k}(1,:) = vpos;
  dat.time{k} = 1/250:1/250:length(vpath{k})/250;
end

cfg=[];
cfg.taper='hanning';
cfg.method = 'mtmfft';
cfg.output='pow';
cfg.foi = 0.01:0.01:10;
cfg.keeptrials = 'yes';
f = ft_freqanalysis(cfg, dat)

for k=1:55
    vpath_sim{k} = simulateVpath(vpath{k},hmmT{k},K);
end

for l=1:10
for k=1:length(vpath)
  vpath_perm{k} = zeros(size(vpath{k}));
  perm_ordering = randperm(K);
  for i=1:K
    vpath_perm{k}(vpath{k}==i)=perm_ordering(i);
  end
end

datperm{l}=[];
datperm{l}.label{1} = 'tmp';
datperm{l}.dimord = '{rpt}_chan_time';
datperm{l}.fsample=250;
for k=1:55
  vpos=nan(length(vpath_perm{k}),1);
for ik=1:12
  vpos(vpath_perm{k}==ik) = disttoplot_manual(ik);
end
  datperm{l}.trial{k}(1,:) = vpos;
  datperm{l}.time{k} = 1/250:1/250:length(vpath{k})/250;
end
fperm{l} = ft_freqanalysis(cfg, datperm{l});
end
  