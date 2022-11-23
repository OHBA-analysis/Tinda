%% Figure 1 Supplement: Plot each figure separately with power and coherence maps


parc=config.parc;
nparcels=config.parc.n_parcels;
local_clim=1;
% cmap = colormap('inferno');
cmap = hotcold;

mni_coords = config.parc.roi_centers;
if whichstudy==4
  alpha = 0.0000001;
else
  alpha = 0.05;
end
if whichstudy==1
  statelabels = {'default mode', 'parietal alpha', 'fronto-temporal', 'visual', 'language', 'default mode', 'sensorimotor', 'parietal alpha', 'dorsal attention', 'R auditory', 'fronto-temporal', 'sensorimotor'};
end

for whichstate =1:K
  
  fig = setup_figure([],2,0.75); axis off
  suptitle(sprintf('RSN state %d - %s', whichstate, statelabels{whichstate}))
  % TINDA
  ax(11) = axes('Position', [0.35 0.74 0.3 0.145]); axis off
  title('TINDA', 'FontSize', 10)
  ax(1) = axes('Position', [0.325 0.5, 0.35, 0.35]);
  cyclicalstateplot_perstate(bestseq,mean_direction,pvals<(alpha/bonf_ncomparisons),find(bestseq==whichstate),false);
  
  ax(9) = axes('Position', [0.05 0.53, 0.25, 0.4]);cla; hold on
  pow = (squeeze(nanmean((psd(:,whichstate,:,:)),4)));
  
  plot(f,(powAvg_freq), '--', 'Color', [.8 .8 .8])
  shadedErrorBar(f,mean(pow,1), std(pow,[],1)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k'},1)
  pow_state_freq{whichstate} = pow;
  set(gca, 'YTick', [])
  set(gca, 'Xtick', nearest(f,30)/4:nearest(f,30)/4:nearest(f,30))
  xlabel('Frequency (Hz)')
  ylabel('PSD')
  [yl(1) yl(2)] = bounds((squash(squeeze(nanmean(nanmean((psd(:,:,nearest(f,2):end,:)),1),4)))));
  ylim([.9 1.1].*yl)
  box off
  
  ax(2) = axes('Position',[0    0.25  0.21 0.21]); % top left
  ax(3) = axes('Position',[0.23  0.25  0.21 0.21]); % top right
  ax(4) = axes('Position',[0.03  0.05 0.21 0.21]);% bottom left
  ax(5) = axes('Position',[0.2   0.05 0.21 0.21]); % bottom right
  
  pow_topo = squeeze(nanmean(nanmean((psd(:,whichstate,:,:)),3),1));
  pow_state_topo{whichstate} = pow_topo;
  toplot = (pow_topo)./(powAvg_topo) - 1;%-mean(net_mean,2);
  if local_clim
    CL = max(abs(toplot(:)))*[-1, 1];%[min(squash(toplot(:,:))) max(squash(toplot(:,:)))];
  else
    CL = max(abs(squash((squeeze(nanmean(nanmean((psd(:,:,:,:)),3),1))))))*[-1, 1];%[min(squash(net_mean(:,:))) max(squash(net_mean(:,:)))]
  end
  psdthresh=CL(1)-.1;%min(net_mean(:))*0.9; % lowest value
  f2 = plot_surface_4way(parc,toplot,0,false,'trilinear',[],psdthresh,CL,ax(2:5));
  
  % coherence topo
  ax(6) = axes('Position',[0.5 0.225  0.24 0.24]);cla
  ax(7) = axes('Position',[0.75  0.225  0.24 0.24]);cla
  ax(8) = axes('Position',[0.625 0.05 0.24 0.24]);cla;
  graph = squeeze(mean(mean(coh(:,whichstate,:,:,:),3),1));
  coh_state_topo{whichstate} = graph;
  [~, ax(6:8), ~] = plot_coh_topo(ax(6:8), mni_coords, graph, cohAvg_topo, [], [], 95);
  
  
  ax(10) = axes('Position', [0.73 0.53, 0.25, 0.4]);cla; hold on
  C = coh(:,:,:,offdiagselect);
  C = (squeeze(nanmean(C(:,whichstate,:,:),4)));
  plot(f, (cohAvg_freq), '--', 'Color', [.8 .8 .8])
  shadedErrorBar(f,mean(C,1), std(C,[],1)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k'},1)
  coh_state_freq{whichstate} = C;
  set(gca, 'YTick', [])
  set(gca, 'Xtick', nearest(f,30)/4:nearest(f,30)/4:nearest(f,30))
  xlabel('Frequency (Hz)')
  ylabel('Coherence')
  [yl(1) yl(2)] = bounds((squash(squeeze(nanmean(nanmean(coh(:,:,nearest(f,2):end,offdiagselect),1),4)))));
  ylim([0.9, 1.1].*yl)
  box off
  
  % colormap for power
  for ii=2:5
    colormap(ax(ii), cmap)
  end
  set_font(10, {'title', 'label'})
  save_figure(fig, [config.figdir, 'figure_supp_tinda_states/', '1supp_tinda_state',int2str(whichstate)],false);
  
  % also save the one with relative x axis
  axes(ax(9))
  cla
  plot(sqrtf, powAvg_freq, '--', 'Color', [.8 .8 .8])
  shadedErrorBar(sqrtf,mean(pow,1), std(pow,[],1)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k'},1)
  set_sqrt_ax(f)
  xlim(sqrtf([1 end]))
  
  axes(ax(10))
  cla
  plot(sqrtf, cohAvg_freq, '--', 'Color', [.8 .8 .8])
  shadedErrorBar(sqrtf,mean(C,1), std(C,[],1)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k'},1)
  set_sqrt_ax(f)
  xlim(sqrtf([1 end]))
  save_figure(fig, [config.figdir,  'figure_supp_tinda_states/', '1supp_tinda_state',int2str(whichstate), '_relative'],false);
  
end

%% Fig 1 supplement: plot as multiple individual state plots:
fig=setup_figure([],2,1);
if whichstudy<4
  cyclicalstateplot_perstate(bestseq,mean_direction,pvals<(0.05/bonf_ncomparisons),[],false,color_scheme);
else
  cyclicalstateplot_perstate(bestseq,mean_direction,pvals<0.0000001*(0.05/bonf_ncomparisons),[],false, color_scheme);
end
save_figure([config.figdir, 'figure_supp_tinda_states/', '1supp_StatePathways']);

% also print for legend:
figure('Position', [440 579 114 219]);
quiver(0,0,1,0,'Color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.8);hold on;
quiver(0,1,1,0,'Color',[0 0 0]+0.8,'LineWidth',2,'MaxHeadSize',0.8);hold on;
axis off;
print([config.figdir, 'figure_supp_tinda_states/', '1supp_StatePathways_legend'], '-dpng')