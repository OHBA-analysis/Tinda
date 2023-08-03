function [fig,ax] = plot_cycle_rt(cfg, bestseq, X, color_scheme)

fig = setup_figure([],2,.9);
Xtmp = zeros(size(X));
Xtmp = X(:,bestseq);
X=Xtmp;
try smoothing=cfg.smoothing; catch,  smoothing = false; end
try cmap=cfg.cmap; catch, cmap=inferno; end
%% polarplot
ax(1) = axes('Position', [.08775 .075 .725 .8]);

if smoothing
  X = imgaussfilt(X, smoothing, 'Padding', 'circular');
end
polarplot3d(flipud(circshift([X],3,2)), 'PolarDirection', 'ccw','AngularRange', [0 2*pi]-pi/12,'TickSpacing', 360,'GridStyle', '-', 'PolarGrid', {0 0},'MeshScale', [1 1]-1./((size(X)./[1 1])),'InterpMethod', 'nearest');

if isfield(cfg, 'clim')
  CL = cfg.clim;
  caxis(CL);
elseif sign(min(X(:))) ~= sign(max(X(:)))
  CL = [-1 1]* max(abs(X(:)));
  caxis(CL);
else
  CL = caxis;
end
view([0,90])
axis off
box off
% colormap((brewermap(256, '')))
colormap(cmap)
%% colorbar
ax(4) = axes('Position', [.65 .175 .3 .6]);
tmp=imagesc(X, CL);
cb = colorbar;
tmp.Visible='off';
box off, axis off
cb.Label.String = cfg.cblabel;
cb.Label.FontSize = 12;
cb.FontSize = 12;
% cb.Limits = [min(X(:)), max(X(:))];

%% cycle plot
ax(2) = axes('Position', [-0.05 0.05 1 .85]);
cyclicalstateplot(bestseq, zeros(12),zeros(12), color_scheme, false)
title(cfg.title)

%% time axis
% extra circle indicating t=0
if cfg.timeaxis(2) == -cfg.timeaxis(1)
  ax(1) = axes('Position', [.08775 .075 .725 .8]);
  r=.475;
  % Create a vectortheta.
  theta=linspace(0,2*pi,200);

  % Generate x-coordinates.
  x=r*cos(theta);

  % Generate y-coordinate.
  y=r*sin(theta);

  % plot the circle.
  plot(x,y,':k', 'LineWidth',2);
  box off, axis off
  xlim([-1 1])
  ylim([-1 1])
end

if cfg.timeaxis(2)>cfg.timeaxis(1)
  ax(3) = axes('Position',[.105 0.475 0.345, 0.001]), plot(cfg.timeaxis(1):cfg.timeaxis(2), 0:0);
else
  ax(3) = axes('Position',[.45 0.475 0.345, 0.001]), plot(cfg.timeaxis(1):cfg.timeaxis(2), 0:0);

end
tmp=ax(3).XTickLabel(1);
ax(3).XTick(1) = 0.95*ax(3).XTick(1);
ax(3).XTickLabel(1) = tmp;
ax(3).FontSize=12;
% ax(3).Color = 'w';
% set(gca,'XColor',[1 1 1]);
xlabel(cfg.timelabel)