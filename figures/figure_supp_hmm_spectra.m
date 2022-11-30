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
      xlim(sqrt([f(1), f(end)]))
    else
      plot(f,Pmean(i,:),'Color',color_scheme{i},'LineStyle',ls{1+(i>6)},'LineWidth',2);
      xlim([f(1), f(end)])
    end
    hold on;
    leglabels{i} = ['State ',int2str(i)];
  end
  
  xlabel('Freq (Hz)'), ylabel('PSD');
  yticks([])
  legend(leglabels, 'Location', 'NorthOutside', 'NumColumns', K/4)
  set_font(10, {'title', 'label'})
  save_figure([config.figdir, 'figure_supp_hmm_spectra/','1supp_PSDperstate', sup]);
end

%% Make a power vs coherence plot

fig = setup_figure([],1,.75); hold on
for k=1:K
  if use_WB_nnmf
    scatter((squeeze(log10(nanmean((psd_wb(:,k,:)),3)))), log10(squeeze(nanmean(coh_wb(:,k,offdiagselect),3))),15, 'MarkerFaceColor', color_scheme{k}, 'MarkerEdgeColor', 'None', 'MarkerFaceAlpha', 0.7);
  else
    scatter((squeeze(log10(nanmean(nanmean((psd(:,k,:,:)),4),3)))), log10(squeeze(nanmean(nanmean((coh(:,k,:,offdiagselect)),4),3))),15, 'MarkerFaceColor', color_scheme{k}, 'MarkerEdgeColor', 'None', 'MarkerFaceAlpha', 0.7);
  end
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
if use_WB_nnmf
  xlim(log10([min(min((squeeze(nanmean((psd_wb),3)))))*0.95, max(max((squeeze(nanmean((psd_wb),3)))))*1.05]))
  ylim(log10([min(min((squeeze(nanmean((coh_wb(:,:,offdiagselect)),3)))))*0.95, max(max((squeeze(nanmean((coh_wb(:,:,offdiagselect)),3)))))*1.05]))
else
  xlim(log10([min(min((squeeze(nanmean(nanmean((psd),4),3)))))*0.95, max(max((squeeze(nanmean(nanmean((psd),4),3)))))*1.05]))
  ylim(log10([min(min((squeeze(nanmean(nanmean((coh(:,:,:,offdiagselect)),4),3)))))*.95, max(max((squeeze(nanmean(nanmean((coh(:,:,:,offdiagselect)),4),3)))))*1.05]))
end
save_figure([config.figdir, 'figure_supp_hmm_spectra/','1supp_PSDvsCoh'] )

%% Do the same for each parcel (averaged over subjects)

fig = setup_figure([],1,.75); hold on
for k=1:K
  if use_WB_nnmf
    scatter((squeeze(log10(nanmean((psd_wb(:,k,:)),1)))), log10(squeeze(nanmean(nanmean(coh_wb(:,k,:,:),1),4))),15, 'MarkerFaceColor', color_scheme{k}, 'MarkerEdgeColor', 'None', 'MarkerFaceAlpha', 0.7);
  else
    scatter((squeeze(log10(nanmean(nanmean((psd(:,k,:,:)),1),3)))), log10(squeeze(nanmean(nanmean(nanmean((coh(:,k,:,:,:)),5),3),1))),15, 'MarkerFaceColor', color_scheme{k}, 'MarkerEdgeColor', 'None', 'MarkerFaceAlpha', 0.7);
  end
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
if use_WB_nnmf
  xlim(log10([min(min((squeeze(nanmean((psd_wb),1)))))*0.95, max(max((squeeze(nanmean((psd_wb),1)))))*1.05]))
  ylim(log10([min(min((squeeze(nanmean(nanmean((coh_wb),4),1)))))*.95, max(max((squeeze(nanmean(nanmean((coh_wb),4),1)))))*1.05]))
else
  xlim(log10([min(min((squeeze(nanmean(nanmean((psd),4),3)))))*0.95, max(max((squeeze(nanmean(nanmean((psd),4),3)))))*1.05]))
  ylim(log10([min(min((squeeze(nanmean(nanmean(nanmean((coh),5),1),3)))))*0.95, max(max((squeeze(nanmean(nanmean(nanmean((coh),5),1),3)))))*1.05]))
end
save_figure([config.figdir, 'figure_supp_hmm_spectra/','1supp_PSDvsCoh_parc'] )


%% Also plot the NNMF if used
if use_WB_nnmf
  fig = setup_figure([],1,1);
  plot(f_orig, wb_comp);
  ylabel('Coherence')
  xlabel('Frequency (Hz)')
  title('2-mode NNMF denoises spectra')
  legend({'mode 1', 'mode 2'})
  save_figure([config.figdir, 'figure_supp_hmm_spectra/','1supp_WB_nnmf'] )
end

