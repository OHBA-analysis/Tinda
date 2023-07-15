addpath('/ohba/pi/mwoolrich/mvanes/scripts/polarplot3d/')
circle_position = zeros(K,2);
for i=1:K
  temp = exp(sqrt(-1)*(i+2)/K*2*pi);
  circle_position(bestseq(i),:) = [real(temp),imag(temp)];
end
circle_position = circle_position(:,1) + sqrt(-1)*circle_position(:,2);
load([config.resultsdir, 'gamma'])
load([config.resultsdir, 'wakehen_events_badseg_adapted.mat']);
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
          trldef{k} = [trldef{k}; k, vis_on, butt_press, rt_tmp, vpath_vis_on, vpath_button];
          vpath_vis_on_all{k} = [vpath_vis_on_all{k}; transpose(vpath_ses{k}(vis_on-125:vis_on+125))];
          vpath_mot_on_all{k} = [vpath_mot_on_all{k}; transpose(vpath_ses{k}(butt_press-250:butt_press))];
          gamma_mot_on_all{k} = [gamma_mot_on_all{k}; squash(double(gamma{k}(butt_press-125:butt_press+125,:)))'];
          gamma_vis_on_all{k} = [gamma_vis_on_all{k}; squash(double(gamma{k}(vis_on-125:vis_on+125,:)))'];
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
save(config.metricfile, "mean_rt_state", '-append')

%% Plot evoked response
meangamma = cellfun(@mean,gamma,'UniformOutput',false);
meangamma = cat(1,meangamma{:});

for k=1:114
  corr_gamma_mot_rt(:,:,k) = reshape(corr(reshape(gamma_mot_on_all{k}, [], 251*12), trldef{k}(:,4)),251,12);
  corr_gamma_vis_rt(:,:,k) = reshape(corr(reshape(gamma_vis_on_all{k}, [], 251*12), trldef{k}(:,4)),251,12);
  gamma_mot_on_all{k} = squeeze(nanmean(gamma_mot_on_all{k}));
  gamma_vis_on_all{k} = squeeze(nanmean(gamma_vis_on_all{k}));
end
gamma_erf_mot = squeeze(nanmean(cat(3,gamma_mot_on_all{:}),3));
gamma_erf_vis = squeeze(nanmean(cat(3,gamma_vis_on_all{:}),3));
%%
fig=setup_figure([],2,1);
subplot(2,2,1), 
plot(-.5:1/250:.5, gamma_erf_vis), colormap(cat(1,color_scheme{:}))
title({'Evoked Gamma', 'timelocked to stimulus onset'})
xlabel('Time (s)'),ylabel('Probability')
for k=1:12
  lg{k} = ['State ', num2str(k)];
end
legend(lg, 'Location', 'NorthWest', 'NumColumns',2)

subplot(2,2,2), 
plot(-1:1/250:0, gamma_erf_mot), colormap(cat(1,color_scheme{:}))
title({'Evoked Gamma', 'timelocked to button press'})
xlabel('Time (s)'), ylabel('Probability')

subplot(2,2,3), 
plot(-.5:1/250:.5, mean(corr_gamma_vis_rt,3)), colormap(cat(1,color_scheme{:}))
title({'Correlation evoked Gamma - RT', 'timelocked to stimulus onset'})
xlabel('Time (s)'), ylabel('R')


subplot(2,2,4), 
plot(-1:1/250:0, mean(corr_gamma_mot_rt,3)), colormap(cat(1,color_scheme{:}))
title({'Correlation evoked Gamma - RT', 'timelocked to button press'})
xlabel('Time (s)'),  ylabel('R')


save_figure([config.figdir, 'figure4_correlations/' 'fig_evoked_gamma'], false)


%% RT circle
X=squeeze(nanmean(mean_rt_state.button.rt(:,bestseq,:)))';
cfg=[];
cfg.timeaxis = mean_rt_state.button.time;
cfg.timelabel = {'Time (s)', 'relative to button press'};
cfg.title = {'Mean RT when in state {\it i}', ' {\it t} seconds before the button press',''};
cfg.cblabel = 'Mean RT (s)';
plot_cycle_rt(cfg, bestseq, X, color_scheme)
save_figure([config.figdir, 'figure4_correlations/' 'fig_rt_cycle_button_press'], false)
%
X=squeeze(nanmean(mean_rt_state.image.rt(:,bestseq,:)))';
cfg=[];
cfg.timeaxis = mean_rt_state.image.time;
cfg.timelabel = {'Time (s)','relative to image onset'};
cfg.title = {'Mean RT when in state {\it i}', ' {\it t} seconds before/after image onset',''};
cfg.cblabel = 'Mean RT (s)';
plot_cycle_rt(cfg, bestseq, X, color_scheme)
save_figure([config.figdir, 'figure4_correlations/' 'fig_rt_cycle_image_onset'], false)

%% Gamma evoked circle


X=gamma_erf_mot(:, bestseq)-mean(meangamma(:,bestseq));
cfg=[];
cfg.timeaxis = [-1 0];
cfg.timelabel = {'Time (s)', 'relative to button press'};
cfg.title = {'Relative state probability',''};
cfg.cblabel = 'Relative probability';
plot_cycle_rt(cfg, bestseq, X, color_scheme)
save_figure([config.figdir, 'figure4_correlations/' 'fig_cycle_gamma_erf_button'], false)
%
X=gamma_erf_vis(:, bestseq)-mean(meangamma(:,bestseq));
cfg=[];
cfg.timeaxis = [-.5 .5];
cfg.timelabel = {'Time (s)','relative to image onset'};
cfg.title = {'Relative state probability',''};
cfg.cblabel = 'Relative probability';
plot_cycle_rt(cfg, bestseq, X, color_scheme)
save_figure([config.figdir, 'figure4_correlations/' 'fig_cycle_gamma_erf_image'], false)

%% Correlation Gamma - RT circle

X=mean(corr_gamma_mot_rt(:,bestseq,:),3);
cfg=[];
cfg.timeaxis = [-1 0];
cfg.timelabel = {'Time (s)', 'relative to button press'};
cfg.title = {'Correlation RT - state probability',''};
cfg.cblabel = 'Correlation';
plot_cycle_rt(cfg, bestseq, X, color_scheme)
save_figure([config.figdir, 'figure4_correlations/' 'fig_cycle_rt_gamma_button'], false)
%
X=mean(corr_gamma_vis_rt(:,bestseq,:),3);
cfg=[];
cfg.timeaxis = [-.5 .5];
cfg.timelabel = {'Time (s)','relative to image onset'};
cfg.title = {'Correlation RT - state probability',''};
cfg.cblabel = 'Correlation';
plot_cycle_rt(cfg, bestseq, X, color_scheme)
save_figure([config.figdir, 'figure4_correlations/' 'fig_cycle_rt_gamma_image'], false)






