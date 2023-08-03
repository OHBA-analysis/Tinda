load([config.resultsdir, 'gamma'])
% reorder gammas
load([config.resultsdir, 'study1matched_state_ordering_MT_nnmf.mat'])
for k=1:114
  gamma{k} = gamma{k}(:, new_state_ordering);
end
load([config.resultsdir, 'wakehen_events_badseg_adapted.mat']);

%%
circle_position = zeros(K,2);
for i=1:K
  temp = exp(sqrt(-1)*(i+2)/K*2*pi);
  circle_position(bestseq(i),:) = [real(temp),imag(temp)];
end
circle_position = circle_position(:,1) + sqrt(-1)*circle_position(:,2);

mean_rt_state.image.time = [-125 125]./config.sample_rate;
mean_rt_state.button.time = [-250 0]./config.sample_rate;

for k=1:length(events)
  k
  trldef{k} = [];
  vpath_vis_on_all{k}=[];
  vpath_mot_on_all{k}=[];
  gamma_mot_on_all{k}=[];
  gamma_vis_on_all{k}=[];
  ix = find(any(events{k}(:,3)==mot_events,2));
  for j=1:length(ix)
    if ix(j)>1
      isvis = any(events{k}(ix(j)-1,3)==vis_events);
      if isvis
        vis_on = double(events{k}(ix(j)-1,1));
        butt_press = double(events{k}(ix(j),1));
        if butt_press+125<length(vpath_ses{k})
          rt_tmp = double(butt_press-vis_on)./config.sample_rate;
          vpath_vis_on = vpath_ses{k}(vis_on);
          vpath_button = vpath_ses{k}(butt_press);
          trldef{k} = [trldef{k}; j, vis_on, butt_press, rt_tmp, vpath_vis_on, vpath_button];
          vpath_vis_on_all{k} = [vpath_vis_on_all{k}; transpose(vpath_ses{k}(vis_on-125:vis_on+125))];
          gamma_vis_on_all{k} = [gamma_vis_on_all{k}; squash(double(gamma{k}(vis_on-125:vis_on+125,:)))'];
          vpath_mot_on_all{k} = [vpath_mot_on_all{k}; transpose(vpath_ses{k}(butt_press-250:butt_press))];
          gamma_mot_on_all{k} = [gamma_mot_on_all{k}; squash(double(gamma{k}(butt_press-250:butt_press,:)))'];
        end
      end
    end
  end
gamma_mot_on_all{k} = reshape(gamma_mot_on_all{k}, [],251, 12);
gamma_vis_on_all{k} = reshape(gamma_vis_on_all{k}, [],251, 12);
  % get average reaction time for each state
  rt = trldef{k}(:,4);
  nperm=100;
  for istate = 1:K
    mean_rt_state.image.rt(k,istate,:) = ((vpath_vis_on_all{k}==istate)'*rt)./sum(vpath_vis_on_all{k}==istate)';
    mean_rt_state.button.rt(k,istate,:) = ((vpath_mot_on_all{k}==istate)'*rt)./sum(vpath_mot_on_all{k}==istate)';
    for iperm=1:nperm
      mean_rt_state.image.rt_perm(k,istate,:,iperm) = ((vpath_vis_on_all{k}==istate)'*rt(randperm(length(rt))))./sum(vpath_vis_on_all{k}==istate)';
      mean_rt_state.button.rt_perm(k,istate,:,iperm) = ((vpath_mot_on_all{k}==istate)'*rt(randperm(length(rt))))./sum(vpath_mot_on_all{k}==istate)';
    end
  end
  % now correlate the reaction times with angle at the moment of vis/mot
  % event
  %{
  for j=1:251
    R_mot(k,j) = circ_corrcl(angle(circle_position(vpath_mot_on_all{k}(:,j))), trldef{k}(:,4));
    R_vis(k,j) = circ_corrcl(angle(circle_position(vpath_vis_on_all{k}(:,j))), trldef{k}(:,4));
    % also for permutations
    nperm=1000;
    for iperm=1:nperm
      perm = randperm(length(trldef{k}(:,4)));
      R_mot_perm(k,j,iperm) = circ_corrcl(angle(circle_position(vpath_mot_on_all{k}(:,j))), trldef{k}(perm,4));
      R_vis_perm(k,j,iperm) = circ_corrcl(angle(circle_position(vpath_vis_on_all{k}(:,j))), trldef{k}(perm,4));
    end
  end
  %}
end

%% Plot evoked response
meangamma = cellfun(@mean,gamma,'UniformOutput',false);
meangamma = cat(1,meangamma{:});

for k=1:114
  for j=1:12
    corr_gamma_mot_rt(:,j,k) = corr(gamma_mot_on_all{k}(:,:,j), trldef{k}(:,4));
  corr_gamma_vis_rt(:,j,k) = corr(gamma_vis_on_all{k}(:,:,j), trldef{k}(:,4));

  for iperm=1:100
    corr_gamma_mot_rt_perm(:,j,k,iperm) = corr(gamma_mot_on_all{k}(:,:,j), trldef{k}(randperm(length(trldef{k})),4));
    corr_gamma_vis_rt_perm(:,j,k,iperm) = corr(gamma_vis_on_all{k}(:,:,j), trldef{k}(randperm(length(trldef{k})),4));
  end
  end
%     gamma_mot_on_all{k} = squeeze(nanmean(gamma_mot_on_all{k}));
%   gamma_vis_on_all{k} = squeeze(nanmean(gamma_vis_on_all{k}));
end
corr_gamma_rt.image.corr = corr_gamma_vis_rt;
corr_gamma_rt.image.time = [-.5 .5];
corr_gamma_rt.button.corr = corr_gamma_mot_rt;
corr_gamma_rt.button.time = [-1 0];
% save(config.metricfile, "mean_rt_state",'corr_gamma_rt', '-append')


gamma_erf_mot = squeeze(nanmean(cat(3,gamma_mot_on_all{:}),3));
gamma_erf_vis = squeeze(nanmean(cat(3,gamma_vis_on_all{:}),3));
%%
tmp = hsv(256);
tmp = tmp(1:17:end,:);
for k=1:12
color_scheme{k} = tmp(k,:);
end
fig=setup_figure([],2,1);
subplot(2,2,1), 
p=plot(-.5:1/250:.5, gamma_erf_vis,'LineWidth',2);
for k=1:12
  p(k).Color = color_scheme{k};
end
colormap(cat(1,color_scheme{:}))
title({'Evoked Gamma', 'timelocked to stimulus onset'})
xlabel('Time (s)'),ylabel('Probability')
for k=1:12
  lg{k} = ['State ', num2str(k)];
end
l=legend(lg, 'Location', 'NorthWest', 'NumColumns',2)
l.Box='off';

subplot(2,2,2), 
p=plot(-1:1/250:0, gamma_erf_mot,'LineWidth',2);
for k=1:12
   p(k).Color = color_scheme{k};
end
title({'Evoked Gamma', 'timelocked to button press'})
xlabel('Time (s)'), ylabel('Probability')

subplot(2,2,3), 
p=plot(-.5:1/250:.5, mean(corr_gamma_vis_rt,3),'LineWidth',2);
% colormap(cat(1,color_scheme{:}));
for k=1:12
   p(k).Color = color_scheme{k};
end
title({'Correlation evoked Gamma - RT', 'timelocked to stimulus onset'})
xlabel('Time (s)'), ylabel('R')


subplot(2,2,4), 
p=plot(-1:1/250:0, mean(corr_gamma_mot_rt,3),'LineWidth',2);
for k=1:12
   p(k).Color = color_scheme{k};
end
title({'Correlation evoked Gamma - RT', 'timelocked to button press'})
xlabel('Time (s)'),  ylabel('R')


save_figure([config.figdir, 'figure4_correlations/' 'fig_evoked_gamma'], false)


%% RT circle
X=squeeze(nanmean(mean_rt_state.button.rt(:,:,:)))';
cfg=[];
cfg.timeaxis = mean_rt_state.button.time;
cfg.timelabel = {'Time (s)', 'relative to button press'};
cfg.title = {'Mean RT when in state {\it i}', ' {\it t} seconds before the button press',''};
cfg.cblabel = 'Mean RT (s)';
plot_cycle_rt(cfg, bestseq, X, color_scheme)
save_figure([config.figdir, 'figure4_correlations/' 'fig_rt_cycle_button_press'], false)
%
X=squeeze(nanmean(mean_rt_state.image.rt(:,:,:)))';
cfg=[];
cfg.timeaxis = mean_rt_state.image.time;
cfg.timelabel = {'Time (s)','relative to image onset'};
cfg.title = {'Mean RT when in state {\it i}', ' {\it t} seconds before/after image onset',''};
cfg.cblabel = 'Mean RT (s)';
plot_cycle_rt(cfg, bestseq, X, color_scheme)
save_figure([config.figdir, 'figure4_correlations/' 'fig_rt_cycle_image_onset'], false)

%% Gamma evoked circle


X=gamma_erf_mot-mean(meangamma);
cfg=[];
cfg.timeaxis = [-1 0];
cfg.timelabel = {'Time (s)', 'relative to button press'};
cfg.title = {'Relative state probability',''};
cfg.cblabel = 'Relative probability';
plot_cycle_rt(cfg, bestseq, X, color_scheme)
save_figure([config.figdir, 'figure4_correlations/' 'fig_cycle_gamma_erf_button'], false)
%
X=gamma_erf_vis-mean(meangamma);
cfg=[];
cfg.timeaxis = [-.5 .5];
cfg.timelabel = {'Time (s)','relative to image onset'};
cfg.title = {'Relative state probability',''};
cfg.cblabel = 'Relative probability';
plot_cycle_rt(cfg, bestseq, X, color_scheme)
save_figure([config.figdir, 'figure4_correlations/' 'fig_cycle_gamma_erf_image'], false)

%% Correlation Gamma - RT circle

X=mean(corr_gamma_mot_rt(:,:,:),3);
cfg=[];
% cfg.clim=[-0.1 0.1];
cfg.smoothing = [5,2];
cfg.cmap = flipud(brewermap(128, 'RdBu'));
cfg.timeaxis = [-1 0];
cfg.timelabel = {'Time (s)', 'relative to button press'};
cfg.title = {''};
cfg.cblabel = 'Correlation';
plot_cycle_rt(cfg, bestseq, X, color_scheme)

save_figure([config.figdir, 'figure4_correlations/' 'fig_cycle_rt_gamma_button'], false)
%
X=mean(corr_gamma_vis_rt(:,:,:),3);
cfg=[];
% cfg.clim=[-0.01 0.01];
cfg.smoothing = [6 2];
cfg.timeaxis = [-.5 .5];
cfg.cmap = flipud(brewermap(128, 'RdBu'));
cfg.timelabel = {'Time (s)','relative to image onset'};
cfg.title = {''};
cfg.cblabel = 'Correlation';
[fig,ax]=plot_cycle_rt(cfg, bestseq, X, color_scheme);
%
save_figure([config.figdir, 'figure4_correlations/' 'fig_cycle_rt_gamma_image'], false)

%%

Nsub=19;
cfgs=[];
cfgs.method = 'montecarlo';
cfgs.statistic='depsamplesT';
cfgs.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfgs.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfgs.tail = 0;
cfgs.correctm = 'cluster';
cfgs.clustertail = 0;

neighbours = [];
for k=1:12
  ix=find(bestseq==k);
  neighbours(k).label = dat.label{k};
  if ix==1
    neighbours(k).neighblabel = dat.label(bestseq([ix+1 12]));
  elseif ix==12
    neighbours(k).neighblabel = dat.label(bestseq([ix-1 1]));
  else
    neighbours(k).neighblabel = dat.label(bestseq([ix-1 ix+1]));
  end
end
cfgs.neighbours = neighbours;
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.numrandomization=1000;
cfgs.parameter='trial';
cfgs.alpha=0.05;
cfgs.clusteralpha = 0.05;

%% button

dat=[];
for k=1:12
  dat.label{k} = ['state', num2str(k)];
end
dat.time=[-1:1/250:0];
dat.dimord = 'rpt_chan_time';
dat.trial = permute(squeeze(mean(reshape(corr_gamma_mot_rt, [251,12,6,19]),3)), [3,2,1]);

datperm=dat;
datperm.trial = permute(squeeze(mean(reshape(mean(corr_gamma_mot_rt_perm,4), [251,12,6,19]),3)), [3,2,1]);
stat=ft_timelockstatistics(cfgs, dat, datperm);
%%

X=stat.stat'.*stat.mask';%mean(corr_gamma_mot_rt(:,:,:),3).*stat.mask';
cfg=[];
% cfg.clim=[-0.1 0.1];
cfg.cmap = flipud(brewermap(128, 'RdBu'));
cfg.timeaxis = [-1 0];
cfg.timelabel = {'Time (s)', 'relative to button press'};
cfg.title = {''};
cfg.cblabel = 't-stat';
plot_cycle_rt(cfg, bestseq, X, color_scheme)
save_figure([config.figdir, 'figure4_correlations/' 'fig_cycle_rt_gamma_button_stat'], false)


%%
dat.trial = permute(squeeze(mean(reshape(corr_gamma_vis_rt, [251,12,6,19]),3)), [3,2,1]);

datperm=dat;
datperm.trial = permute(squeeze(mean(reshape(mean(corr_gamma_vis_rt_perm,4), [251,12,6,19]),3)), [3,2,1]);
cfgs.time=[-.5:1/250:.5];
stat=ft_timelockstatistics(cfgs, dat, datperm);


X=stat.stat'.*stat.mask';%mean(corr_gamma_mot_rt(:,:,:),3).*stat.mask';
cfg=[];
% cfg.clim=[-0.1 0.1];
cfg.cmap = flipud(brewermap(128, 'RdBu'));
cfg.timeaxis = [-.5 .5];
cfg.timelabel = {'Time (s)', 'relative to button press'};
cfg.title = {''};
cfg.cblabel = 't-stat';
plot_cycle_rt(cfg, bestseq, X, color_scheme)
save_figure([config.figdir, 'figure4_correlations/' 'fig_cycle_rt_gamma_image_stat'], false)
