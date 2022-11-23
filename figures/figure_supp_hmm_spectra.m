%% plot average state spectra vs freq (averaged over all parcels):

Pmean = squeeze(mean(mean(psd,1),4));

% relorabs = 'rel'; % 'rel' or 'abs'
for relorabs = {'abs', 'rel'}
  fig = setup_figure([],1,1);
  if strcmp(relorabs, 'abs')
    sup = '';
  elseif strcmp(relorabs, 'rel')
    sup = '_relative';
  end
  ls = {'-','--'};
  for i=1:hmm.K
    if strcmp(relorabs, 'rel')
      plot(sqrtf,Pmean(i,:),'Color',color_scheme{i},'LineStyle',ls{1+(i>6)},'LineWidth',2);
      set_sqrt_ax(f)
      xlim(sqrt([f(1), 30]))
    else
      plot(f,Pmean(i,:),'Color',color_scheme{i},'LineStyle',ls{1+(i>6)},'LineWidth',2);
      xlim([f(1), 30])
    end
    hold on;
    leglabels{i} = ['State ',int2str(i)];
  end
  
  xlabel('Freq (Hz)'), ylabel('PSD');
  yticks([])
  legend(leglabels, 'Location', 'NorthOutside', 'NumColumns', K/4)
  set_font(10, {'title', 'label'})
  save_figure([config.figdir,'1supp_PSDperstate', sup]);
end

%% Make a power vs coherence plot

fig = setup_figure([],1,.75); hold on
for k=1:K
  scatter((squeeze(log10(nanmean(nanmean((psd(:,k,:,:)),4),3)))), log10(squeeze(nanmean(nanmean((coh(:,k,:,offdiagselect)),4),3))),15, 'MarkerFaceColor', color_scheme{k}, 'MarkerEdgeColor', 'None', 'MarkerFaceAlpha', 0.7);
  l{k} = sprintf('State %d', k);
end
% axis off
box off
yticks([])
ylabel('Coherence')
xticks([])
xlabel('PSD')
legend(l, 'Location', 'EastOutside', 'NumColumns', 1)
box off
axis square
set_font(10, {'label', 'title'})
xlim(log10([min(min((squeeze(nanmean(nanmean((psd),4),3)))))*0.95, max(max((squeeze(nanmean(nanmean((psd),4),3)))))*1.05]))
ylim(log10([min(min((squeeze(nanmean(nanmean((coh(:,:,:,offdiagselect)),4),3)))))*1.05, max(max((squeeze(nanmean(nanmean((coh(:,:,:,offdiagselect)),4),3)))))*0.95]))

save_figure([config.figdir,'1supp_PSDvsCoh'] )
