%% Figure 2 supplement:  analyse by quintiles
percentiles = 0:20:100;
fig=setup_figure([],2,2);
if whichstudy==4
    alpha_thresh = 0.0000001;
else
    alpha_thresh = 0.05;
end
clear ax
for ip=1:length(percentiles)-1
    [FO_p,FO_pvals_p,t_intervals_p, FO_stat_p] = computeLongTermAsymmetry(vpath,hmmT,K,percentiles(ip:ip+1));
    hmm_1stlevel.control.quintile.assym_permtest{ip} = FO_stat_p;
    hmm_1stlevel.control.quintile.assym_permtest{ip}.pvals = FO_pvals_p;
    hmm_1stlevel.control.quintile.t_intervals{ip} = t_intervals_p;
    a=[];
    for i=1:K
        for j=1:K
            [a.h(i,j), a.pvals(i,j), a.ci(i,j,:), a.stat(i,j)] = ttest(squeeze(FO_p(i,j,1,:)), squeeze(FO_p(i,j,2,:)));
        end
    end
    hmm_1stlevel.control.quintile.assym_ttest{ip} = a;
    t_int(ip,:) = round(mean(mean(cellfun(@mean, t_intervals_p),2)));
    mean_direction_p = squeeze(mean(FO_p(:,:,1,:)-FO_p(:,:,2,:),4));
    FO_assym_p = squeeze((FO_p(:,:,1,:)-FO_p(:,:,2,:))./mean(FO_p,3));
    %   bestseq_p = optimiseSequentialPattern(FO_p);
    %   bestseq_p = bestseq_p{1};
    %   angleplot_p = circle_angles(bestseq_p);
    bestseq_p = bestseq;
    angleplot_p = circle_angles(bestseq_p);
    
    rotational_momentum_p = compute_rotational_momentum(angleplot_p, FO_assym_p); squeeze(imag(nansum(nansum(angleplot_p.*FO_assym_p))));
    tmp(:,ip) = rotational_momentum_p;
    
    ax(ip,1) = axes('Position', [-.05, 1.02-0.2*ip, 0.45, 0.13]);
    
    cyclicalstatesubplot(bestseq_p,mean_direction_p,a.pvals<hmm_1stlevel.assym_ttest.alpha_thresh);
    title({sprintf('Rotational momentum = %0.3f', mean(rotational_momentum_p, 'omitnan')/hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum),''})
    
    
    ax(ip,2) = axes('Position', [0.45, 1.03-0.196*ip, 0.5, 0.15]);
    ITmerged = cellfun(@mean,t_intervals_p);ITmerged = ITmerged./config.sample_rate;
    tmp2(ip) = median(mean(ITmerged,2));
    distributionPlot(sqrt(ITmerged),'showMM',2,'color',{color_scheme{1:size(FO,2)}});
    set(ax(ip,2),'YLim',[0 1.1*max(sqrt(ITmerged(:)))])
    if ip==5
        xlabel('RSN-State')
    end
    grid on; ylabel('Interval Times (s)')
    for ii=1:length(ax(ip,2).YTickLabel)
        ax(ip,2).YTickLabel{ii} = num2str(str2num(ax(ip,2).YTickLabel{ii})^2);
    end
    hmm_1stlevel.control.quintile.FO_intervals(:,:,:,ip) = squeeze([FO_p(:,:,1,:) - FO_p(:,:,2,:)]./mean(FO_p,3));
    hmm_1stlevel.control.quintile.bestseq{ip} = bestseq_p;
    hmm_1stlevel.control.quintile.t_intervals{ip} = t_intervals_p;
end
hmm_1stlevel.control.quintile.rotational_momentum = tmp;
set_font(10, {'title', 'label'})

save_figure([config.figdir, 'figure_supp_tinda_quintiles/','2supp_Cyclicalpatterns_percentiled']);


%% Make figure with cycle plots and scatter/boxplots underneath
% fig=setup_figure([], 2, 0.4);
% clear t_int
% for k=1:5
%     ax(k) = axes('Position', [0.035+0.2*(k-1) 0.5, 0.14 0.5]);
%     FO_p = hmm_1stlevel.control.quintile.FO_intervals(:,:,:,k);
%     mean_direction_p = mean(FO_p,3);
%     a = hmm_1stlevel.control.quintile.assym_ttest{k};
%     cyclicalstatesubplot(bestseq,mean_direction_p,a.pvals<hmm_1stlevel.assym_ttest.alpha_thresh);
%     rotational_momentum_p = hmm_1stlevel.control.quintile.rotational_momentum(:,k);
%     t_int(k,:,:) = cellfun(@mean, hmm_1stlevel.control.quintile.t_intervals{k});
% end
% ax(6) = axes('Position', [0.05, 0.05, .9 0.4]);
% 
% A1 = nan(config.nSj,14);
% A2=A1;
% A1(:, 1:3:end) = mean(t_int,3)';
% A2(:, 2:3:end) = hmm_1stlevel.control.quintile.rotational_momentum./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum;
% yyaxis left
% boxplot_with_scatter(A1, clr{1}, 0.5)
% yticks([0 500 1000]), yticklabels({'0', '0.5', '1'})
% 
% 
% yyaxis right
% boxplot_with_scatter(A2, clr{2}, 0.5);
% hold on, h=hline(0),
% box off,
% ylim([-.12 .22])
% xticks(1.5:3:14),  xlabel('Quintile')
% gcf()
% text(.0, .25, 'Mean IT', 'Color', clr{1})
% text(12, .25, 'Rotational Momentum', 'Color', clr{2})
% save_figure([config.figdir, 'figure_supp_tinda_quintiles/','2supp_cycle_quintiled']);

%%
studylabel = {'MEGUK', 'MEGUK', 'HCP', 'Cam-CAN', 'HCP', 'Cam-CAN', 'WakeHen'};

fig=setup_figure([], 1, 3);
ax(8) = axes('Position', [.5 .97 .001, .001]);
box off
axis off
t=title(studylabel{whichstudy});
t.FontSize = 14;

clear t_int
for k=1:5
    ax(k) = axes('Position', [-0.02 .98-(k)*0.19, 0.7 0.15]);
    FO_p = hmm_1stlevel.control.quintile.FO_intervals(:,:,:,k);
    mean_direction_p = mean(FO_p,3);
    a = hmm_1stlevel.control.quintile.assym_ttest{k};
    cyclicalstatesubplot(bestseq,mean_direction_p,a.pvals<hmm_1stlevel.assym_ttest.alpha_thresh);
    rotational_momentum_p = hmm_1stlevel.control.quintile.rotational_momentum(:,k);
    t_int(k,:,:) = cellfun(@mean, hmm_1stlevel.control.quintile.t_intervals{k});
    t = text(-1.3, -0.45, ['Quintile ', num2str(k)]);
    t.Rotation=90;
    t.FontWeight = 'bold';
end


q = hmm_1stlevel.control.quintile.rotational_momentum./hmm_1stlevel.cycle_metrics.max_theoretical_rotational_momentum;
q2 = mean(t_int,3)';

ax(6) = axes('Position', [0.65, 0.03, 0.3 .85]);
ax(6).Color = 'none';
bar(5:-1:1, mean(q), 0.25, 'EdgeColor',clr{2}, 'FaceColor', clr{2}) 

box off
xticks([])
xlim([.8 5.2])
yticks([0 0.05 0.1])

ylabel('Mean M (a.u.)', 'Color', clr{2})
set(gca, 'color', 'none');
view([90, -90])

ax(7) = axes('Position', [0.65, 0.085, 0.3 .85]);
ax(7).Color = 'none';
bar(5:-1:1, mean(q2), 0.25, 'EdgeColor',clr{1}, 'FaceColor', clr{1})
ax(7).YAxisLocation = 'right';
xlim([0.8 5.2])
yticks([100, 500, 1000])

xticks([])
ylabel('Mean IT (ms)', 'Color', clr{1})
view([90, -90])
box off
set(gca, 'color', 'none');

save_figure([config.figdir, 'figure_supp_tinda_quintiles/','2supp_cycle_quintiled']);



%% Make an extra plot with median interval time vs rotational momentum - with smaller sliding window
[FO_p,pvals_p,t_intervals_p] = computeLongTermAsymmetry(vpath,hmmT,K,percentiles(ip:ip+1));

percentiles = 5:5:95;
clear rotational_momentum_p
for ip=1:length(percentiles)
    [FO_p,pvals_p,t_intervals_p] = computeLongTermAsymmetry(vpath,hmmT,K,[percentiles(ip)-4 percentiles(ip)+5]);
    FO_assym_p = squeeze((FO_p(:,:,1,:)-FO_p(:,:,2,:))./mean(FO_p,3));
    ITmean(ip) = mean(mean(cellfun(@mean,t_intervals_p)));ITmean(ip) = ITmean(ip)./config.sample_rate;
    rotational_momentum_p(ip,:) = compute_rotational_momentum(angleplot, FO_assym_p);
    hmm_1stlevel.control.ntile.FO_intervals(:,:,:,ip) = squeeze([FO_p(:,:,1,:) - FO_p(:,:,2,:)]./mean(FO_p,3));
end
hmm_1stlevel.control.ntile.rotational_momentum = rotational_momentum_p;

fig = setup_figure([],2,0.6);
subplot(1,2,1),
shadedErrorBar(percentiles,nanmean(rotational_momentum_p,2), nanstd(rotational_momentum_p,[],2)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k'},1)
xlabel('Percentile'), ylabel('Rotational Momentum')
yticks([0])
subplot(1,2,2)
shadedErrorBar((ITmean),nanmean(rotational_momentum_p,2), nanstd(rotational_momentum_p,[],2)./sqrt(config.nSj), {'LineWidth', 2, 'Color', 'k'},1)
yticks([0])
xlabel('Mean IT (s)'), ylabel('Rotational momentum')
save_figure([config.figdir, 'figure_supp_tinda_quintiles/', '2supp_Cyclicalpatterns_rotationalmomentum_percentiled']);
close