%% Compare the different TINDA measures
clr = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]};

q = hmm_1stlevel.cycle_metrics;
% positive only measures first
tmp = [q.rotational_momentum./q.max_theoretical_rotational_momentum, q.FO_assym_subject_fit, q.TIDA, q.circularity_subject];
ttl = {{'Rotational', 'momentum'}, {'FO asym', 'fit'}, {'TIDA', ''}, {'Circularity', ''}};
l = size(tmp,2);

fig = setup_figure([],2,0.5);
for k=1:l
  ax(k) = axes('Position', [0.05+(k-1)*(.6/l) 0.55 (0.3/l) 0.4],'Color','none','XColor','none'); hold on
  scatter(1+ones(size(tmp(:,k))).*((rand(size(tmp(:,k)))-0.5)/2),tmp(:,k),'filled', 'MarkerFaceAlpha',0.6)
  h = boxplot(tmp(:,k), 'Colors', 'k', 'Width', .75, 'whisker', 1000) % put this on topend
  set(h, {'linew'}, {2})
  title(ttl{k})
  xticks([])
  box off
end

ax(5) = axes('Position', [0.05 0.05 0.525 0.4],'Color','none','XColor','none'); hold on
tmp2 = (tmp./std(tmp,[],1))';
tmp2(3:4,:) = tmp2(3:4,:) - min(tmp2(3:4, :)')';
plot(tmp2, '-o', 'color', 0.8*[1,1,1,0.3])
plot(mean(tmp2,2), '-', 'color', 'k', 'LineWidth', 2)
plot(mean(tmp2,2), '.', 'color', 'k', 'MarkerSize', 30)

xlim([0.66 4.33])
xticks([])
yticks(0)
title('Cycle measures over subjects')

c = flipud(brewermap(64, 'RdBu'));
r = corr(tmp);
r(find(triu(r)))=0;
ax(6) = axes('Position', [0.65 0.12 0.33 0.76]);
imagesc(1:4, 1:4, r, [-1, 1]); colormap(c); colorbar; box off
[R, C] = ndgrid(1:4, 1:4);
R = R(:); C = C(:) - 1/4;
%rows are Y values, columns are X values !
vals = r(:);
vals(vals==0)=nan;
text(C, R, string(round(vals,2)), 'color', 'w')
xlim([.5 3.5])
ylim([1.5 4.5])
xticks(1:3)
yticks(2:4)
% xticklabels(1:3, 'HorizontalAlignment', 'center')
set(gca(),'XTickLabel',{sprintf('  Rot\\newlinemom') sprintf('FO asym\\newline       fit') sprintf('TIDA\\newline ')});
set(gca(),'YTickLabel',{sprintf('  FO\\newlineasym\\newline   fit') sprintf('TIDA\\newline ') sprintf('Circ\\newline ') });
title({'Correlation between' 'cycle measures'})
set_font(10, {'title', 'label'})

save_figure([config.figdir, 'figure_supp_tinda_metrics/','2supp_cycle_metrics'])

%% Do the same for the per - state measures


clr = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]};

q = hmm_1stlevel.cycle_metrics;
% positive only measures first
tmp = cat(3,q.rotational_momentum_perstate./q.max_theoretical_rotational_momentum, q.TIDA_perstate);
ttl = {{'Rotational momentum'}, {'TIDA'}};
l = size(tmp,3);


fig = setup_figure([],2,0.6);
for ik=1:l
  %   ax = axes('Position', [0.05 0.05+(ik-1)*ik])
  %   subplot(2,6,(ik-1)*5+[ik:ik+3]), hold on
  %   subplot(2,2,2*(ik-1)+1)
  ax(2*(ik-1)+1) = axes('Position', [0.1 0.1+0.5*(2-ik) 0.4 0.35]);
  hold on
  for k=1:K
    scatter(1+ones(size(tmp(:,k,ik))).*((k-1)+(rand(size(tmp(:,k,ik)))-0.5)/2),tmp(:,k,ik),'filled', 'MarkerFaceAlpha',0.6, 'MarkerFaceColor',color_scheme{k})
  end
  h = boxplot(tmp(:,:,ik), 'Colors', 'k', 'Width', .75, 'whisker', 1000); % put this on topend
  set(h, {'linew'}, {2})
  title(ttl{ik})
  xlabel('State')
  
  %   subplot(2,2,2*ik)
  ax(2*ik) = axes('Position', [0.575 0.1+0.5*(2-ik) 0.4 0.35], 'YDir', 'reverse');
  hold on
  r = corr(tmp(:,:,ik));
  r(find(triu(r)))=0;
  c2 = c(29:end,:);
  imagesc(1:12, 1:12, r, [0, 1]);
  colormap(c2); colorbar; box off
  [R, C] = ndgrid(1:12, 1:12);
  R = R(:); C = C(:) - 1/3;
  %rows are Y values, columns are X values !
  vals = r(:);
  vals(vals==0)=nan;
  text(C, R, string(round(vals,1)), 'color', 'k', 'FontSize', 6)
  xlim([.5 11.5])
  ylim([1.5 12.5])
  xticks(1:11)
  yticks(2:12)
  title({'Correlation between states'})
  xlabel('State')
  ylabel('State')
end
save_figure([config.figdir, 'figure_supp_tinda_metrics/','2supp_cycle_metrics_perstate'])


%% Compare Cycle measures with simulations
q0 = hmm_1stlevel.cycle_metrics;
q1 = hmm_1stlevel.simulation{1}.cycle_metrics;
q2 = hmm_1stlevel.simulation_average.cycle_metrics;
tmp = [];
tmp(:,:,1) = [q0.rotational_momentum./q0.max_theoretical_rotational_momentum, q1.rotational_momentum./q0.max_theoretical_rotational_momentum, q2.rotational_momentum./q0.max_theoretical_rotational_momentum];
tmp(:,:,2) = [q0.FO_assym_subject_fit, q1.FO_assym_subject_fit, q2.FO_assym_subject_fit];
tmp(:,:,3) = [q0.TIDA, q1.TIDA, q2.TIDA];
tmp(:,:,4) = [q0.circularity_subject, q1.circularity_subject, q2.circularity_subject];

ttl = {{'Rotational momentum'}, {'FO asym fit'}, {'TIDA'}, {'Circularity'}};
measures = {'rotational_momentum', 'FO_assym_subject_fit', 'TIDA', 'circularity_subject'};%, 'TIDA_perstate', 'rotational_momentum_perstate';
yl = [-0.2 0.15; -0.4 1.1; 0 0.25; 0 1];

fig = setup_figure([],2,.4);
for ik = 1:4
  ax(ik) = axes('Position', [0.05 + (ik-1)*0.25, 0.1 0.18 0.8]);
  
  hold on
  for k=1:3
    scatter(1+ones(size(tmp(:,k,ik))).*((k-1)+(rand(size(tmp(:,k,ik)))-0.5)/2),tmp(:,k,ik),'filled', 'MarkerFaceAlpha',0.6);%, 'MarkerFaceColor',color_scheme{k})
  end
  h = boxplot(tmp(:,:,ik), 'Colors', 'k', 'Width', .75, 'whisker', 1000); % put this on topend
  set(h, {'linew'}, {2})
  title(ttl{ik})
  sigstar({[1,2], [1,3]}, [hmm_1stlevel.metric_vs_sim.(measures{ik}).prob, hmm_1stlevel.metric_vs_sim_avg.(measures{ik}).prob])
  xticks([1:3])
  xticklabels({sprintf('observed') sprintf('1 sim') sprintf('100 sim')});
  xtickangle(45)
  
  ylim(yl(ik,:))
  box off
end
save_figure([config.figdir, 'figure_supp_tinda_metrics/','2supp_cycle_metrics_vs_sim'])

%% Compare per state cycle measures with simulations

q0 = hmm_1stlevel.cycle_metrics;
q1 = hmm_1stlevel.simulation{1}.cycle_metrics;
q2 = hmm_1stlevel.simulation_average.cycle_metrics;
tmp = [];
tmp(:,:,:,1) = cat(3, q0.rotational_momentum_perstate./q0.max_theoretical_rotational_momentum_perstate, q1.rotational_momentum_perstate./q0.max_theoretical_rotational_momentum_perstate, q2.rotational_momentum_perstate./q0.max_theoretical_rotational_momentum_perstate);
tmp(:,:,:,2) = cat(3, q0.TIDA_perstate, q1.TIDA_perstate, q2.TIDA_perstate);


ttl = {{'Rotational momentum'}, {'TIDA'}};
measures = {'rotational_momentum_perstate', 'TIDA_perstate'};%, 'TIDA_perstate', 'rotational_momentum_perstate';
if whichstudy==1 || whichstudy==2
  yl = [-0.3 0.35; -0.05 .45];
elseif whichstudy==3
  yl = [-0.3 0.35; -0.05 .4];
elseif whichstudy==4
  yl = [-0.3 0.35; -0.05 .4];
end


for ik = 1:length(measures)
  fig = setup_figure([],2,.6);
  for whichstate=1:12
    if whichstate<=K/2
      ax(whichstate) = axes('Position', [0.05 + (whichstate-1)*0.96/6, 1.05-1*0.5 0.65/6+0.02 0.4]);
    else
      ax(whichstate) = axes('Position', [0.05 + (whichstate-6-1)*0.96/6, 1.05-2*0.5 0.65/6+0.02 0.4]);
    end
    hold on
    
    k=1; scatter(1+ones(size(tmp(:,whichstate, k,ik))).*((k-1)+(rand(size(tmp(:,whichstate, k,ik)))-0.5)/2),tmp(:,whichstate,k,ik),'filled', 'MarkerFaceAlpha',0.3, 'MarkerFaceColor',color_scheme{whichstate})
    k=2; scatter(1+ones(size(tmp(:,whichstate, k,ik))).*((k-1)+(rand(size(tmp(:,whichstate, k,ik)))-0.5)/2),tmp(:,whichstate,k,ik),'d', 'MarkerEdgeAlpha',0.3, 'MarkerEdgeColor',color_scheme{whichstate})
    k=3; scatter(1+ones(size(tmp(:,whichstate, k,ik))).*((k-1)+(rand(size(tmp(:,whichstate, k,ik)))-0.5)/2),tmp(:,whichstate,k,ik),'filled','d', 'MarkerFaceAlpha',0.3, 'MarkerFaceColor',color_scheme{whichstate})
    
    h = boxplot(squeeze(tmp(:,whichstate,:,ik)), 'Colors', 'k', 'Width', .75, 'whisker', 1000); % put this on topend
    set(h, {'linew'}, {2})
    sigstar({[1,2], [1,3]}, [hmm_1stlevel.metric_vs_sim.(measures{ik}).prob(whichstate), hmm_1stlevel.metric_vs_sim_avg.(measures{ik}).prob(whichstate)])

    xticks([1:3])
    xticklabels({sprintf('observed') sprintf('1 sim') sprintf('100 sim')});
    xtickangle(45)
    if whichstate~=1 && whichstate~=K/2+1
      yticks([])
    end
    ylim(yl(ik,:))
    title(sprintf('State %d', whichstate))
    box off
  end
  suptitle(ttl{ik})
end
save_figure([config.figdir, 'figure_supp_tinda_metrics/','2supp_cycle_metrics_vs_sim'])
