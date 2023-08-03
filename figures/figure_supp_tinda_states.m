%% Figure 1 Supplement: Plot each figure separately with power and coherence maps

parc=config.parc;
hotcold = cmap_hotcold;

nparcels=config.parc.n_parcels;
local_clim=1;
cmap = hotcold;

mni_coords = config.parc.roi_centers;
alpha = 0.05;

if whichstudy==1
    th = 95;
else
    th = 98;
end

fig = setup_figure([],2,1.5); axis off
cl{1} = [0 0.4470 0.7410];
cl{2} = [0.8500 0.3250 0.0980];
w=0;q=1;
for k =1:K
    if k>K/2
        w=0.5; % used for placing plots in second column
        q=1.02;
    end
    
    % Power spectrum
    ax(k, 1) = axes('Position', [q*w+0.325, .925-2.1*((k-1)-w*K)/(K+1), .15, .75*1/K]); hold on
    pow = pow_state_freq{k};
    shadedErrorBar(sqrtf,mean(pow,1), std(pow,[],1)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k', 'LineStyle', '-'},1)
    plot(sqrtf, powAvg_freq, '--', 'Color', .8*[1 1 1])
    set_sqrt_ax(f)
    xlim(sqrtf([1 end]))
    set(gca, 'YTick', [])
    xticklabels([])
    ylabel('PSD')
    [yl(1) yl(2)] = bounds((squash(squeeze(nanmean(nanmean((psd(:,:,nearest(f,2):end,:)),1),4)))));
    ylim([.9 1.1].*yl)
    box off
    
    % Coherence spectrum
    ax(k, 2) = axes('Position', [q*w+0.325, .8575-2.1*((k-1)-w*K)/(K+1), .15, 0.75*1/K]); hold on
    C = coh_state_freq{k};
    shadedErrorBar(sqrtf,mean(C,1), std(C,[],1)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k', 'LineStyle', '-'},1)
    plot(sqrtf, cohAvg_freq, '--', 'Color', .8*[1 1 1])
    set_sqrt_ax(f)
    xlim(sqrtf([1 end]))
    set(gca, 'YTick', [])
    xlabel('Frequency (Hz)')
    ylabel('Coherence')
    [yl(1) yl(2)] = bounds((squash(squeeze(nanmean(nanmean(coh(:,:,nearest(f,2):end,offdiagselect),1),4)))));
    ylim([0.9, 1.1].*yl)
    box off
    
    % power map
    ax(k, 3) = axes('Position', [q*w+0.01, .92-2.1*((k-1)-w*K)/(K+1), 1.4*1/K, 1/K]);
    ax(k, 4) = axes('Position', [q*w+0.178, .92-2.1*((k-1)-w*K)/(K+1), 1.4*1/K, 1/K]);
    
    pow_topo =   pow_state_topo{k};
    toplot = (pow_topo)./(powAvg_topo) - 1;
    CL = max(abs(toplot(:)))*[-1, 1];

    psdthresh=CL(1)-.1;
    f2 = plot_surface_4way(parc,toplot,0,false,'trilinear',[],psdthresh,CL,ax(k, 3:4));

    % color bar
    ax(k,10) = axes('Position', [q*w+.116, .912-2.1*((k-1)-w*K)/(K+1), .07, .05]), 
    tmp = imagesc(CL);
    tmp.Visible = 'off'; box off, axis off;
    cb = colorbar('southoutside');
    
    cb.Limits = CL*100;
    cb.Ticks = [CL(1)*100 0 CL(2)*100];
    cb.TickLabels = {num2str(ceil(cb.Ticks(1))), num2str(0), num2str(floor(cb.Ticks(3)))};
    cb.Ruler.TickLabelRotation=0;
    ax(k,10).CLim = 100*CL;

    % coherence map
    ax(k, 6) = axes('Position', [q*w+0.0175, .815-2.1*((k-1)-w*K)/(K+1), 1.2*.8*1/K, 1.4*1/K]); cla, box off, axis off
    ax(k, 7) = axes('Position', [q*w+0.205, .815-2.1*((k-1)-w*K)/(K+1), 1.2*.8*1/K, 1.4*1/K]); cla, box off, axis off
    ax(k, 8) = axes('Position', [q*w+0.115, .83-2.1*((k-1)-w*K)/(K+1), 1.1*.8*1/K, 1.1*1/K]); cla, box off, axis off
    
    graph = coh_state_topo{k};
    [~, ax(6:8), ~] = plot_coh_topo(ax(k,6:8), mni_coords, graph, cohAvg_topo, [1 2], [], th);
    ax(k, 9) = axes('Position', [q*w+0.112, .91-2.1*((k-1)-w*K)/(K+1), 0.075, .75*1/K]); cla, axis off, box off
    title(sprintf('State %d', k))

end
%%
% colormap for power
for k=1:K
    for ii=[3:4 10]
        colormap(ax(k,ii), cmap)
    end

end

save_figure([config.figdir,  'figure_supp_tinda_states/', '1supp_all_state_descriptions', '_relative'],[],false);
close