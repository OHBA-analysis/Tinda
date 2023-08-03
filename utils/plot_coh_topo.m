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
  sphereCols = repmat([0 0 0], nROIs, 1);%repmat([30 144 255]/255, nROIs, 1);
end
if ~exist('perc', 'var') || isempty(perc)
  perc = 95;
end

graph = coh_state - coh_mean;
colorLims = [];

viewZ = {[270,0],[-270,0],[0,90]};

for iplot=1:length(ax)
  axes(ax(iplot))
  osl_braingraph(graph, colorLims, repmat(0.5,nROIs,1), [0 1], mni_coords, [], perc, sphereCols, edgeLims);
  view(ax(iplot),viewZ{iplot})
  colorbar('hide')
end