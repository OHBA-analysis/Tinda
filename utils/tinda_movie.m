function tinda_movie(bestseq, mean_direction, sigpoints, f, psd, coh, outname)
% this function creates a movie from the output of the Tinda analysis and
% the output of loadMTspect (from Higgins 2020, Neuron). The following is
% expected:
% bestseq: 1xN (Tinda sequence)
% mean_direction: NxN (Tinda mean direction)
% sigpoints: NxN bool (the to be plotted Tinda connections)
% f: 1xnfreq frequency bins
% psd: nsub x nstate x nfreq x nparc x nparc OR the output from
% teh_spectral_nnmf containing at least nnmf_psd_* "specs" and "maps"
% coh: nsub x nstate x nfreq x nparc x nparc OR the output from
% teh_spectral_nnmf containing at least nnmf_coh_* "specs" and "maps"
% outname: str



nstates = length(bestseq);
fsel=30;
disttoplot_manual = zeros(nstates,2);
for i=1:nstates
  temp = exp(sqrt(-1)*(i+2)/nstates*2*pi);
  disttoplot_manual(bestseq(i),:) = [real(temp),imag(temp)];
end

doreplay=0;
color_scheme = set1_cols();
cmap = colormap('inferno');
numframes=120;%480;
v=VideoWriter(outname);
v.FrameRate = 7.5;
position = [0 0 1120 840];
open(v)
nrow = 6;
ncol = 6;

% some coherence topo set up
parcelFile = fullfile(osldir,'parcellations',    'fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
parc=parcellation(parcelFile);
spatialRes = 8;
spatialMap = nii.quickread(parcelFile, spatialRes);
nparcels = parc.n_parcels;
mni_coords = find_ROI_centres(spatialMap, spatialRes, 0, osldir);
viewZ = {[270,0],[-270,0],[0,90]};

% reshape PSD
if isstruct(coh) && isstruct(psd)
  psd_freq(1,:,:) = nanmean(reshape(transpose(psd.nnmf_psd_specs)*reshape(permute(psd.nnmf_psd_maps, [2,1,3]),2, []), [], nstates, nparcels),3)';
  psd_freq(2,:,:) = psd_freq(1,:,:); % to keep the right dimensions
  psd_topo = squeeze(psd.nnmf_psd_maps(:,1,:))';
  % coh_freq can't be recreated from the info present (for some reason?)
  % use the "raw" coh instead
  coh_freq = repmat(coh.coh_freq, [2,1,1]);
%   coh_freq(1,:,:) = nanmean(nanmean(reshape(transpose(coh.nnmf_coh_specs)*reshape(permute(coh.nnmf_coh_maps, [2,1,3,4]),2, []), [], nstates ,nparcels, nparcels),4),3)';
  coh_topo = squeeze(coh.nnmf_coh_maps(:,1,:,:));
  psd_freq = psd_freq(:,:,1:nearest(f,30));
  coh_freq = coh_freq(:,:,1:nearest(f,30));
else
  diagselect = find(eye(nparcels));
  notdiagselect = find(~eye(nparcels));
  psd = psd(:,:,1:nearest(f,fsel),diagselect); % take the diagonal of nparc x nparc
  psd_topo = zeros(nparcels,nstates);
  for k = 1:size(psd,2)
    psd_topo(:,k) = squeeze(nanmean(nanmean(psd(:,k,:,:),3),1))';
  end
  psd_freq = nanmean(psd,4);
  
  % reshape coh
  coh = coh(:,:,1:nearest(f,fsel),:,:);
  coh_topo = squeeze(nanmean(nanmean(coh,3), 1));
  coh_freq = squeeze(nanmean(coh(:,:,:,notdiagselect),4));
end

f = f(1:nearest(f,30));


for ivid=[numframes/4:-1:1, numframes:-1:numframes/4]
  %%%%%%%%%
  % CLOCK %
  %%%%%%%%%
  fig=figure('Position',position);

  ax(5)=subplot(nrow,ncol,[3,4,9,10]);
  cyclicalstateplot(bestseq,mean_direction,sigpoints,[],0)
  %plot vector simulated:
  hold on;
  temp = exp(sqrt(-1)*(ivid)/numframes*2*pi);
  plot([0 real(temp)],[0,imag(temp)],'r','LineWidth',2);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CREATE MOVING AVERAGE WEIGHTS %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %     xpdf = mvnpdf(disttoplot_manual,[real(temp),imag(temp)],0.5*eye(2));
%   q = abs((disttoplot_manual(:,1) + sqrt(-1)*disttoplot_manual(:,2)) - temp);
  q = abs((disttoplot_manual) - temp);
  [val, idx]=sort(abs(q), 'ascend');
  xpdf = zeros(nstates,1);
  if val(1)<0.001
    xpdf(idx(1))=1;
  else
    xpdf(idx(1:2)) = 1./val(1:2);
  end
  xpdf = xpdf./sum(xpdf);
  if doreplay
    %%%%%%%%%%%%%%%
    % REPLAY PLOT %
    %%%%%%%%%%%%%%%
    ax(6)=subplot(nrow,ncol, [5,11]);
    bar(sum(replay_counts*xpdf)), ylim([0 1200]), axis off, title('Replay counts'),
    set(gca,'fontsize',16);
  end
  
  %%%%%%%%%%%%
  % PSD PLOT %
  %%%%%%%%%%%%
  psd_temp = squeeze(sum(psd_freq.*repmat(xpdf',[size(psd_freq,1),1,size(psd_freq,3)]),2));
  subplot(nrow,ncol,[15,16,21,22])%11);
  shadedErrorBar(f,mean(psd_temp,1),std(psd_temp,[],1)./sqrt(size(psd_temp,1)),{'-k', 'LineWidth', 2});
  q1=mean(psd_freq) + std(psd_freq);
  yl=[min(q1(:)), max(q1(:))];
  ylim(yl);
  xlim([1,30]);
%   xlabel('Freq (Hz)')
  ylabel('PSD')

  %%%%%%%%%%%%%%%%%%
  % COHERENCE PLOT %
  %%%%%%%%%%%%%%%%%%
  coh_temp = squeeze(sum(coh_freq.*repmat(xpdf',[size(coh_freq,1),1,size(coh_freq,3)]),2));
  subplot(nrow,ncol,[27,28,33,34])%12);
  shadedErrorBar(f,mean(coh_temp,1),std(coh_temp,[],1)./sqrt(size(coh_temp,1)),{'-k', 'LineWidth', 2});
  q1=mean(coh_freq) + std(coh_freq);
  yl=[min(q1(:)), max(q1(:))];
  ylim(yl);
  xlim([1,30]);
  xlabel('Freq (Hz)')
  ylabel('Coherence')
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  % COHERENCE TOPOGRAPHY %
  %%%%%%%%%%%%%%%%%%%%%%%%
%   graph = abs(reshape(xpdf'*reshape(squeeze(nnmfWB_res.nnmf_coh_maps(:,1,:,:)), K, []), nparcels, nparcels));
  graph = reshape(squeeze(xpdf'*reshape(coh_topo, length(xpdf), [])), nparcels, nparcels);
  graph(find(eye(nparcels))) = 0;
  tmp=squash(triu(graph));
  inds2=find(tmp>1e-10);
  data=tmp(inds2);
  S2=[];
  S2.data=squash(data);
  S2.do_fischer_xform=0;
  S2.pvalue_th=0.05/(nparcels.^2);
  S2.do_plots=0;
  graph_ggm=teh_graph_gmm_fit(S2);
  
  th=graph_ggm.normalised_th;
  graph=graph_ggm.data';
  if th<1.96 % less than 2 stds from mean
    graph(graph<th)=NaN;
    graphmat=nan(nparcels, nparcels);
    graphmat(inds2)=graph;
    graph=graphmat;
  else
    % few sparse connections, do not plot:
    graph = nan(nparcels);
  end
  if all(isnan(graph(:)))
    graph(1,1)=1;
  end
  
  nROIs = size(graph,2);
  colorLims = [th th+1];
  sphereCols = repmat([30 144 255]/255, nROIs, 1);
  edgeLims = [2 4];
  
%   figure('Color', 'w','Position',[547 100 577 453]);
  ax(3) = subplot(nrow,ncol,[25,26,31,32]);%axes('Position',[0 0.5 0.5 0.5]);
  ax(4) = subplot(nrow,ncol,[29,30,35,36]);%axes('Position',[0.55 0.5 0.5 0.5]);
%   ax(5) = %axes('Position',[0.27 0.1 0.5 0.5]);
  
  % and plot 3 way brain graphs:
  for iplot=1:2%3
    axes(ax(2+iplot))
    osl_braingraph(graph, colorLims, repmat(0.5,nROIs,1), [0 1], mni_coords, [], 0, sphereCols, edgeLims);
    view(ax(2+iplot),viewZ{iplot})
    colorbar('hide')
  end
    
  %%%%%%%%%%%%%%%%%%%%
  % POWER TOPOGRAPHY %
  %%%%%%%%%%%%%%%%%%%%
  toplot = psd_topo*xpdf;%-mean(net_mean,2);
  if 1 % scale colorlim per power map
    CL = max(abs(squash(toplot(:,:)))) * [0 1];
  else % scale color lim over all power maps
    CL = max(abs(squash(psd_topo(:,:)))) * [0 1];
  end
  if 0 % threshold the power map
    psdthresh = prctile(abs(toplot),50);
    toplot(abs(toplot)<psdthresh) = NaN;
  end

  ax(1)=subplot(nrow,ncol,[13,14,19,20]); axis off
  ax(2)=subplot(nrow, ncol, [17,18,23,24]); axis off

  plot_surface(parc,toplot,0,false,[], ax(1:2), CL);
  colormap(cmap)
  
  %%%%%%%%%%%%%%%
  % WRITE VIDEO %
  %%%%%%%%%%%%%%%
  set(fig, 'Position', position);
  gfig=getframe(fig);
  writeVideo(v,gfig)
  close(gcf)
end
close(v);