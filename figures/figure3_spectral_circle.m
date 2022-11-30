%% Figure 3: Spectral information in the circle plot
local_clim=true;
diffmode = {'rel'}; % whether to plot psd/coh relative to average over states
use_sqrt_f = [true];%[true, false] ; % whether to plot PSD/coh with sqrt of f axis
statecolor = [true, false]; % whether to plot PSD/coh in the state color
parc = config.parc;
mni_coords = parc.template_coordinates;
clear yl tmp*
for whichtopo = 1:2%3
  do_pow_or_coh = {'pow','coh'};% 'both'
  do_pow_or_coh = do_pow_or_coh{whichtopo};
  for sc = statecolor
    if sc
      sup2 = '_statecolor';
    else
      sup2 = '_samecolor';
    end
    for ii = use_sqrt_f
      if ii
        sup1 = '_sqrtax';
      else
        sup1 = '_linax';
      end
      
      pos = zeros(12,2);
      for i=1:12
        temp = exp(sqrt(-1)*(i+2)/12*2*pi);
        pos(i,:) = [real(temp),imag(temp)];
      end
      % manually shift some positions
      pos = (pos/4 + 0.25)*2;
      pos(:,1) = pos(:,1)*0.88;
      pos(pos(:,2)>0.55,2) = pos(pos(:,2)>0.55,2) - 0.05;
      pos(pos(:,2)<0.4,2) = pos(pos(:,2)<0.4,2) + 0.05;
      pos([2 6],1) = pos([2 6],1) - [0.075;0.075];
      pos([8 12],1) = pos([8 12],1) + [0.075;0.075];
      
      pos(bestseq,:) = pos;
      fig = setup_figure([],2,1);
      ax(13,1) = axes('Position', [0.295 0.19 0.38 .6]); hold on
      axes(ax(13,1))
      cyclicalstateplot(bestseq,hmm_1stlevel.cycle_metrics.mean_direction, zeros(12), color_scheme, false)
      
      for k=1:12
        ax(k,2) = axes('Position', [0.125 0.025 0 0]+[0.85 0.85 1 1].*[pos(k,1), pos(k,2), 0.1 0.1]); hold on
        yyaxis left
        yyaxis right
        if sc
          cl{1} = color_scheme{k};
          cl{2} = cl{1};
          linestyle = ':';
          ax(k,2).YAxis(1).Color = 'k';
          ax(k,2).YAxis(2).Color = 'k';
        else
          cl{1} = [0 0.4470 0.7410];
          cl{2} = [0.8500 0.3250 0.0980];
          linestyle = '-';
        end
        
        yyaxis left
        tmpP = mean(pow_state_freq{k});
        if strcmp(diffmode{1}, 'rel')
          tmpP = (tmpP./powAvg_freq) - 1;
          [ylt(1) ylt(2)] = bounds(squash((squeeze(mean(cat(3,pow_state_freq{:}))))./powAvg_freq'-1));
          [ylt(3) ylt(4)] = bounds(squash((squeeze(mean(cat(3,coh_state_freq{:}))))./cohAvg_freq'-1));
          [yl(1) yl(2)] = bounds(ylt);
        else
          [yl(1) yl(2)] = bounds(squash(squeeze(nanmean(nanmean((psd(:,:,nearest(f,3):end,:)),1),4))));
        end
        
        if ii
          plot(sqrtf,tmpP, 'LineWidth',2, 'Color', cl{1});
          set_sqrt_ax(f)
        else
          plot(f,tmpP, 'LineWidth',2, 'Color', cl{1});
        end
        hline(0, ':k')
        
        % ylim([0.95,1.05].*yl)
        ylim([1.1,1.1].*yl)
        
        yticks([])
        if sc
          ylabel(['PSD ' char(8212)]);
        else
          ylabel('PSD');
        end
        yyaxis right
        tmpC = mean(coh_state_freq{k});
        if strcmp(diffmode{1}, 'rel')
          tmpC = (tmpC./cohAvg_freq)-1;
        else
          [yl(1) yl(2)] = bounds(squash(squeeze(nanmean(nanmean((coh(:,:,nearest(f,3):end,offdiagselect)),1),4))));
        end
        
        if ii
          plot(sqrtf,tmpC, 'LineWidth',2, 'Color', cl{2}, 'LineStyle', linestyle);
          xlim(sqrt([f(1) 30]))
          set_sqrt_ax(f)
        else
          plot(f,tmpC, 'LineWidth',2, 'Color', cl{2}, 'LineStyle', linestyle);
          xlim([f(1) 30])
        end
        ylim([1.1 1.1].*yl)
        hline(0, ':k')
        
        if sc
          ylabel(['Coh ' '---']);
        else
          ylabel('Coh');
        end
        yticks([])
        box off
        
        
        % pow topo
        toplot = pow_state_topo{k}./powAvg_topo-1;
        coords_left = mni_coords(:,1)<0;
        %       thresh = prctile(abs(toplot),0);
        % project the right side to the left side (all non thresholded values are shown on the left projection)
        tmp = [toplot(1:2:end), toplot(2:2:end)];
        [~, ix] = max(abs([toplot(1:2:end), toplot(2:2:end)]),[],2, 'omitnan');
        for ii=1:length(ix)
          tmp2(ii) = tmp(ii,ix(ii));
        end
        toplot(1:2:end) = tmp2;
        toplot(2:2:end) = tmp2;
        %       toplot(toplot<thresh)=NaN;
        if local_clim
          CL = 1.1*max(abs(pow_state_topo{k}./powAvg_topo - 1));
          CL = [-CL CL];
        else
          CL = 1.1*max(abs(squash(cat(2,pow_state_topo{:})./powAvg_topo)-1));
          CL = [-CL CL];
        end
        
        
        % Coh topo
        graph = coh_state_topo{k};
        
        if strcmp(do_pow_or_coh, 'both')
          ax(k,1) = axes('Position', [0 0.065 0 0]+[0.85 0.85 1 1].*[pos(k,1), pos(k,2), 0.1 0.1]);
          plot_surface_4way(parc,toplot,1,true,'trilinear', [],CL(1)-0.1,CL, ax(k,1))

          ax(k,3) = axes('Position', [0.01 -0.015 0 0]+[0.85 0.85 1 1].*[pos(k,1), pos(k,2), 0.08 0.08]);
          [~, ax(k,3), ~] = plot_coh_topo(ax(k,3), mni_coords, graph, cohAvg_topo, [0 3], [], 95);  axis off      
        elseif strcmp(do_pow_or_coh, 'pow')
          ax(k,1) = axes('Position', [0 0.025 0 0]+[0.85 0.85 1 1].*[pos(k,1), pos(k,2), 0.1 0.1]);
          plot_surface_4way(parc,toplot,1,true,'trilinear', [],CL(1)-0.1,CL, ax(k,1))
        elseif  strcmp(do_pow_or_coh, 'coh')
          ax(k,1) = axes('Position', [0.01 0.025 0 0]+[0.85 0.85 1 1].*[pos(k,1), pos(k,2), 0.08 0.08]);
          [~, ax(k,1), ~] = plot_coh_topo(ax(k,1), mni_coords, graph, cohAvg_topo, [0 3], [], 95);axis off
        end
        colormap(hotcold)
      end
      ax(11,1) = axes('Position', [0.365, 0.38, 0.25, 0.25]); hold on
      clear l
      for k=1:K
        scatter(log10(squeeze(nanmean(nanmean((psd(:,k,:,:)),4),3))), log10(squeeze(nanmean(nanmean(coh(:,k,:,offdiagselect),4),3))),15, 'MarkerFaceColor', color_scheme{k}, 'MarkerEdgeColor', 'None', 'MarkerFaceAlpha', 0.7);
        l{k} = sprintf('State %d', k);
      end
      % axis off
      box off
      yticks([])
      ylabel('Coherence')
      xticks([])
      xlabel('PSD')
      box off
      axis square
      xlim(log10([min(min((squeeze(nanmean(nanmean((psd),4),3)))))*0.95, max(max((squeeze(nanmean(nanmean((psd),4),3)))))*1.05]))
      ylim(log10([min(min((squeeze(nanmean(nanmean(coh(:,:,:,offdiagselect),4),3)))))*1.05, max(max((squeeze(nanmean(nanmean(coh(:,:,:,offdiagselect),4),3)))))*0.95]))
      set_font(10, {'label', 'title'})
      save_figure([config.figdir, 'figure3_spectral_circle/','3_Spectral_circle_', do_pow_or_coh, sup1, sup2], false)
    end
  end
end