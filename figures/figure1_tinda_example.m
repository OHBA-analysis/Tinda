%% Figure 1: Plot TINDA example
if whichstudy==1
  iSj=5;
  if strcmp(config.reordering_states, 'coherence')
    whichstate=1;
  else
    whichstate=2;
  end
  fprintf(['\n Subject: ',int2str(iSj), ', State: ' int2str(whichstate)]);
  
  % load raw data
  D = spm_eeg_load(replace(hmm.data_files{iSj}, '/Users/chiggins/data/YunzheData/Replaydata4Cam/WooliePipeline/spm/', '/ohba/pi/mwoolrich/datasets/ReplayData/Neuron2020Analysis/'));
  
  % get vpath
  Gamma_subj = Gamma(hmm.subj_inds==iSj,:);
  [~,vpath_subj] = max(Gamma_subj,[],2);
  
  
  % parcellation
  parc=config.parc;
  mni_coords = config.parc.roi_centers;
  nparcels = length(parc.labels);
  
  % Figure setup
  fig = setup_figure([],2,0.6); clear ax
  cmap = colormap('inferno');
  local_clim = true; % for the power maps, create a clim based on that state's range
  
  %%%%%%%%%%%%%%%%%%%
  % DATA TRACE PLOT %
  %%%%%%%%%%%%%%%%%%%
  ax(1) = axes('Position', [.1 .55 .25 .4]); hold on
  
  tlength=720;
  if (whichstate==1 && strcmp(config.reordering_states, 'replay')) || (whichstate==2 && strcmp(config.reordering_states, 'coherence'))
    tstart=34300;
    % this is the interval we'll plot TINDA for
    t1 = 106;
    t2 = 668;
  elseif (whichstate==2 && strcmp(config.reordering_states, 'replay')) || (whichstate==1 && strcmp(config.reordering_states, 'coherence'))
    tstart=30300;
    t1 = 76;
    t2 = 697;
  end
  t_segment = tstart + [1:tlength];
  thalf = mean([t1, t2]);
  t_segment1 = tstart + [t1:tlength/2];
  t_segment2 = tstart + [tlength/2+1:t2];
  
  q = vpath_subj(t_segment)==whichstate;
  s1 = [find(diff(q)==1)+1 find(diff(q)==-1)];
  mn = min(min(D(1:8,t_segment)));
  mx = max(max(D(1:8,t_segment)));
  
  %   plot_Gamma(Gamma_subj(t_segment,:), t_segment, true,false,[]), colororder(cat(1,color_scheme{:}))
  %   a = gca;
  %   for k=1:length(a.Children)
  %     a.Children(k).FaceAlpha=0.2;
  %     a.Children(k).LineWidth=0.001;
  %   end
  %   box off
  %   ylim([.02 .99])
  
  % plot raw trace
  plot((1/250):(1/250):(length(t_segment)/250), D(1:8,t_segment(1:tlength))')
  % annotate visits to whichstate
  for jj=1:size(s1,1)
    fill([s1(jj,1) s1(jj,1) s1(jj,2) s1(jj,2)]./250, 1.3*[mn mx mx mn], color_scheme{whichstate}, 'FaceAlpha', .5, 'EdgeColor', 'none');
  end
  set(ax(1),'xcolor','none')
  ax(1).XAxis.Label.Color=[0 0 0];
  ax(1).XAxis.Label.Visible='on';
  ax(1).YAxis.Label.Color=[0 0 0];
  set(ax(1),'YTick',[]);
  axis off
  text(ax(1).Position(1)-0.3, ax(1).Position(2)+0.05, {'Resting', 'state', 'data'}, 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10)
  
  %%%%%%%%%%%%%%
  % VPATH PLOT %
  %%%%%%%%%%%%%%
  ax(2) = axes('Position', [0.1 0.1 0.25 0.4]);
  hold on
  
  % annotate intervals
  cb = [256,193,1; 201, 204, 231]/256; % interval annotation colors
  fill([t1-1 t1-1 thalf thalf]./250, [0 13 13 0], cb(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none') % interval T1
  fill([thalf thalf t2 t2]./250, [0 13 13 0], cb(2,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none')% interval T2
  
  yl = [0 13];
  % plot vpath
  for k=1:12
    myline = NaN(length(t_segment),1);
    myline(vpath_subj(t_segment)==k)=(13-k);
    plot((1/250):(1/250):(length(t_segment)/250),myline,'LineWidth',20,'Color',color_scheme{k});
    yaxislabels{13-k} = ['State ',int2str(k)];
  end
  set(ax(2),'YTick',[1:12]);
  yaxislabels{13-whichstate} = ['\color{red} ' yaxislabels{13-whichstate}];
  set(ax(2),'YTickLabels',yaxislabels);
  grid off
  set(ax(2),'xcolor','none')
  ax(2).XAxis.Label.Color=[0 0 0];
  ax(2).XAxis.Label.Visible='on';
  xlabel('Time');
  ax(2).YAxis.Label.Color=[0 0 0];
  ylim(yl)
  title('HMM state timecourse', 'HorizontalAlignment', 'center')
  %   title('FO', 'HorizontalAlignment', 'center')
  %   title('\boldmath$FO^T$', 'Interpreter', 'Latex', 'FontSize', 10)
  
  text(0.28, -.06, sprintf('%s', sym("T_1")), 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold')
  text(0.72, -.06, sprintf('%s', sym("T_2")), 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold')
  
  
  %%%%%%%%%%%%%%%%
  % T1 BAR GRAPH %
  %%%%%%%%%%%%%%%%
  ax(3) = axes('Position', [0.375 0.1 0.075 0.4]);
  hold on
  
  
  for ii=1:12
    T1(ii) = sum(vpath_subj(t_segment1)==ii);
    h=barh(ii, T1(ii));
    set(h,'FaceColor', color_scheme{ii})
  end
  set(gca, 'Ydir', 'reverse')
  axis off
  ylim(yl)
  %   title({'$$\ \ \ \ \ \ \ \forall_{k\neq{1}}\sum_{T_{1}|T_{2}}Z_{t}==k$$'},'Interpreter', 'Latex', 'FontSize', 10, 'HorizontalAlignment', 'center')
  title('\boldmath$FO^{T_1}$', 'Interpreter', 'Latex', 'FontSize', 10, 'FontWeight', 'bold')
  
  %%%%%%%%%%%%%%%%
  % T2 BAR GRAPH %
  %%%%%%%%%%%%%%%%
  ax(4) = axes('Position', [0.475 0.1 0.075 0.4]);
  hold on
  for ii=1:12
    T2(ii) = sum(vpath_subj(t_segment2)==ii);
    h=barh(ii, T2(ii));
    set(h,'FaceColor', color_scheme{ii})
  end
  set(gca, 'Ydir', 'reverse')
  axis off
  ylim(yl)
  title('\boldmath$FO^{T_2}$', 'Interpreter', 'Latex', 'FontSize', 10, 'FontWeight', 'bold')
  
  %%%%%%%%%%%%%%%%%%%%%
  % DISTRIBUTION PLOT %
  %%%%%%%%%%%%%%%%%%%%%
  sigpoints = hmm_1stlevel.FO_pvals<(alpha_thresh/bonf_ncomparisons);
  for ii=1:12
    ax(4+ii) = axes('Position', [0.575 0.121+(12-ii)*0.031 0.1 0.022]);
    hold on
  end
  
  clear d
  cb = [256,193,1; 201, 204, 231]/256;
  for ii=1:12
    axes(ax(4+ii))
    if ii==1
      %         title({'$$\forall_{k\neq{1}} \sum_{i=1}^{N}Q_{i,T_{1}}-Q_{i,T_{2}}$$'},'Interpreter', 'Latex', 'FontSize', 10, 'HorizontalAlignment', 'center')
      title('FO asym', 'FontSize', 10)
    end
    if any(setdiff(1:12, whichstate)==ii)
      d{1} = squeeze(hmm_1stlevel.FO_intervals(whichstate,ii,1,:));
      d{2} = squeeze(hmm_1stlevel.FO_intervals(whichstate,ii,2,:));
      % exagerate significant differences
      if sigpoints(ii,whichstate)
        if hmm_1stlevel.cycle_metrics.mean_direction(ii,whichstate)<0
          d{1} = d{1}+0.5*mean(d{2});
        elseif hmm_1stlevel.cycle_metrics.mean_direction(ii,whichstate)>0
          d{2} = d{2}+0.5*mean(d{2});
        end
      end
      h2 = raincloud_plot(d{2}, 'box_on', 0, 'color', cb(2,:), 'alpha', 0.5);
      h1 = raincloud_plot(d{1}, 'box_on', 0, 'color', cb(1,:), 'alpha', 0.5);
      y2 = get(ax(4+ii), 'YLim');
      ylim([0 y2(2)*1.3])
    end
    axis off
    set(gca, 'YTick', []);
    set(gca, 'XTick', []);
  end
  
  %%%%%%%%%%%%%%%%%%%%%
  % TINDA CIRCLE PLOT %
  %%%%%%%%%%%%%%%%%%%%%
  ax(25) = axes('Position', [0.70 0.3500 0.3 0.145]); axis off
  title('TINDA', 'FontSize', 10)
  ax(17) = axes('Position', [0.70 0.0500 0.3 0.4]);
  cyclicalstateplot_perstate(bestseq,hmm_1stlevel.cycle_metrics.mean_direction,sigpoints,find(bestseq==whichstate),false);
  % seq=[12:-1:1];
%     cyclicalstateplot_perstate(seq,mean_direction,sigpoints,find(seq==whichstate),false);
  
  
  %%%%%%%%%%%
  % PSD MAP %
  %%%%%%%%%%%
  ax(18) = axes('Position',[0.34        0.8 0.175 0.2] ); % top left
  ax(19) = axes('Position',[0.34+0.185  0.8  0.175 0.2] ); % top right
  ax(20) = axes('Position',[0.34+0.025  0.6  0.175 0.2] ); % bottom left
  ax(21) = axes('Position',[0.34+0.16   0.6  0.175 0.2] ); % bottom right
  ax(26) = axes('Position', [0.37 0.80 0.3 0.15]); axis off
  title('PSD', 'FontSize', 10)
    pow_topo = pow_state_topo{whichstate};

  toplot = (pow_topo)./(powAvg_topo) - 1;%-mean(net_mean,2);
  if local_clim
    CL = max(abs(toplot(:)))*[-1, 1];%[min(squash(toplot(:,:))) max(squash(toplot(:,:)))];
  else
    CL = max(abs(squash((squeeze(nanmean(nanmean((psd(:,:,:,:)),3),1))))))*[-1, 1];%[min(squash(net_mean(:,:))) max(squash(net_mean(:,:)))]
  end
  psdthresh=CL(1)-.1;%min(net_mean(:))*0.9; % lowest value
  f2 = plot_surface_4way(parc,toplot,0,false,'trilinear',[],psdthresh,CL,ax(18:21));
  
  %%%%%%%%%%%%%%%%%
  % COHERENCE MAP %
  %%%%%%%%%%%%%%%%%
  ax(22) = axes('Position',[0.63+0.02 0.765 0.24 0.24] ); % left
  ax(23) = axes('Position',[0.63+0.18  0.765 0.24 0.24] ); % right
  ax(24) = axes('Position',[.74 0.6 0.22 0.22]);% bottom
  
  graph = coh_state_topo{whichstate};
  [~, ax(22:24), ~] = plot_coh_topo(ax(22:24), mni_coords, graph, cohAvg_topo, [], [], 95);

  ax(27) = axes('Position', [0.7 0.80 0.3 0.15]); axis off
  title('Coh', 'FontSize', 10, 'HorizontalAlignment', 'center')
  
  
  % change colormap for power
  for ii=18:21
    colormap(ax(ii), hotcold)
  end
  
  %%%%%%%%%%%%
  % SAVE FIG %
  %%%%%%%%%%%%%
  set_font(10, {'label', 'title'})
  save_figure(fig, [config.figdir, 'figure1_tinda_example/', '/1_tinda_example_state',int2str(whichstate)]);
  %   close
end