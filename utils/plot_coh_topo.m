function [fig, ax, graph] = plot_coh_topo(ax, mni_coords, coh_state, coh_mean, edgeLims, sphereCols, perc)

nparcels = size(coh_state);
nROIs = nparcels(1);
if exist('ax', 'var') && ~isempty(ax)
  fig = gcf;
else
  fig = figure;
  for k=1:3
    ax(k) = subplot(1,3,k);
  end
end
if ~exist('coh_mean', 'var') || isempty(coh_mean)
  coh_mean = zeros(nROIs);
end
if ~exist('edgeLims', 'var') || isempty(edgeLims)
  edgeLims = [2 4];
end
if ~exist('sphereCols', 'var') || isempty(sphereCols)
  sphereCols = repmat([30 144 255]/255, nROIs, 1);
end
if ~exist('perc', 'var') || isempty(perc)
  perc = 95;
end

tmp=squash(triu(coh_state-coh_mean));
inds2=find(tmp>1e-10);
if isempty(inds2)
  graph=nan(nparcels);
  graph(1,1)=1;
  colorLims = [0, 1];
else
  data= tmp(inds2);
  
  th = percentile(data,perc);
  graph = nan(nparcels);
  graph(inds2(find(data>=th))) = data(find(data>=th));
  if max(graph(:))==th
    colorLims = [0 th];
%   elseif max(graph(:),[],'omitnan')>2*th
%     colorLims = [th 0.8*max(graph(:), [], 'omitnan')];
  else
    colorLims = [th max(graph(:),[],'omitnan')];
  end
end


%   % if there's more than 50 connections, cut off
%   [val, ix] = sort(graph(:), 'descend', 'MissingPlacement', 'last');
%   if ~isnan(val(51))
%     graph(ix(51:end))=NaN;
%   end
%   S2=[];
%   S2.data=squash(data);
%   S2.do_fischer_xform=0;
%   S2.pvalue_th=0.05/length(S2.data);%(nparcels.^2);
%   S2.do_plots=0;
% %     a=teh_graph_gmm_fit(S2);
% % sum(a.data>a.normalised_th)
%   graph_ggm=teh_graph_gmm_fit(S2);
%
%   th=graph_ggm.normalised_th;
% %   th = percentile(graph_ggm.data, 95);
%   graph=graph_ggm.data';
%
%   if th<1.96 % less than 2 stds from mean
%     graph(graph<th)=NaN;
%     graphmat=nan(nparcels, nparcels);
%     graphmat(inds2)=graph;
%     graph=graphmat;
%   else
%     % few sparse connections, do not plot:
%     graph = nan(nparcels);
%   end
%
%   if all(isnan(graph(:)))
%     graph(1,1)=1;
%   end
%       [~, tmp] = sort(graph(:), 'descend', 'MissingPlacement', 'last');
%       graph(tmp(round(nparcels.^2/50):end)) = NaN;

%   colorLims = [th th+1];


viewZ = {[270,0],[-270,0],[0,90]};

for iplot=1:length(ax)
  axes(ax(iplot))
  osl_braingraph(graph, colorLims, repmat(0.5,nROIs,1), [0 1], mni_coords, [], 0, sphereCols, edgeLims);
  view(ax(iplot),viewZ{iplot})
  colorbar('hide')
end