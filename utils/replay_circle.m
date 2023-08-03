load([config.resultsdir, 'replay_evoked_gamma_perc1.mat'])
% 
load([config.resultsdir, 'replay_ordering.mat'])

G_evoked = reorder_data(replay_ordering,new_state_ordering,G_evoked);

figure; subplot(1,2,1), p=plot(mean(G_evoked,3)); 
cs = colorscheme(1);
for k=1:12
  p(k).Color = cs{k};
end

subplot(1,2,2), p=plot(mean(G2,3)); 
cs = colorscheme(1);
for k=1:12
  p(k).Color = cs{k};
end

%%
X=round(100*mean(G_evoked,3));
cfg=[];
% cfg.clim=[-0.1 0.1];
cfg.smoothing = [5,.5];
cfg.cmap = flipud(brewermap(128, 'RdBu'));
cfg.timeaxis = [-.5 1];
cfg.timelabel = {'Time (s)', 'relative to replay onset'};
cfg.title = {'Memory replay is embedded in cycle phase',''};
cfg.cblabel = 'State probability increase (%)';
[fig,ax]=plot_cycle_rt(cfg, replay.bestseq.modelhmm, X, color_scheme);
set(ax(2), 'FontSize',14);
save_figure([config.figdir, 'figure4_correlations/replay_evoked_gamma_circle'],false)

%% get behavior
d=dir('/ohba/pi/mwoolrich/datasets/ReplayData/ReplayData4Cam/BehavData/*/*Replay*');
for k=1:21
  q=load([d(k).folder,'/' d(k).name]);
  acc1(k,1) = mean(q.data.Subdata.StructLearn.Results(:,4));
end

d=dir('/ohba/pi/mwoolrich/datasets/ReplayData/StrLearn_MEGexp/BehavData/*/*Str*');
for k=1:22
  q=load([d(k).folder,'/' d(k).name]);
  acc2(k,1) = mean(q.data.Subdata.StructLearn.Results(:,4));
end

% correlations
for k=1:12
    R1(k,:,1) = corr(acc1, squeeze(G_evoked(:,k,1:2:42))', 'type', 'Spearman');
    R1(k,:,2) = corr(acc1, squeeze(G_evoked(:,k,2:2:42))', 'type', 'Spearman');
    R2(k,:,1) = corr(acc2, squeeze(G_evoked(:,k,43:2:end))', 'type', 'Spearman');
    R2(k,:,2) = corr(acc2, squeeze(G_evoked(:,k,44:2:end))', 'type', 'Spearman');
    R3(k,:) = corr([acc1;acc1;acc2;acc2], squeeze(G_evoked(:,k,:))', 'type', 'Spearman');

    for iperm=1:1000
      rnd1 = randperm(21);
      rnd2 = randperm(22);
    R1_perm(k,:,1,iperm) = corr(acc1(rnd1), squeeze(G_evoked(:,k,1:2:42))', 'type', 'Spearman');
    R1_perm(k,:,2,iperm) = corr(acc1(rnd1), squeeze(G_evoked(:,k,2:2:42))', 'type', 'Spearman');
    R2_perm(k,:,1,iperm) = corr(acc2(rnd2), squeeze(G_evoked(:,k,43:2:end))', 'type', 'Spearman');
    R2_perm(k,:,2,iperm) = corr(acc2(rnd2), squeeze(G_evoked(:,k,44:2:end))', 'type', 'Spearman');
    R3_perm(k,:,iperm) = corr([acc1(rnd1);acc1(rnd1);acc2(rnd2);acc2(rnd2)], squeeze(G_evoked(:,k,:))', 'type', 'Spearman');
    
    end
end



%%
X=mean(R2,3)';%-mean(mean(R2_perm,3),4)';%transpose((mean(R1,3)+mean(R2,3))./2);

cfg=[];
% cfg.clim=[-0.1 0.1];
% cfg.smoothing = [5,1];
cfg.cmap = flipud(brewermap(128,'RdBu'));
cfg.timeaxis = [-.5 .5];
cfg.timelabel = {'Time (s)', 'relative to replay onset'};
cfg.title = {'Correlation RT - state probability',''};
cfg.cblabel = 'Correlation';
plot_cycle_rt(cfg, replay.bestseq.modelhmm, X, color_scheme)




