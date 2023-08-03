function [fig,ax] = plot_cycle_rt(cfg, bestseq, X, color_scheme)

fig = setup_figure([],2,.9);
Xtmp = zeros(size(X));
Xtmp = X(:,bestseq);
X=Xtmp;
try smoothing=cfg.smoothing; catch,  smoothing = [5,1]; end
%% polarplot
ax(1) = axes('Position', [.08775 .075 .725 .8]);
X = imgaussfilt(X, smoothing, 'Padding', 'circular');
polarplot3d(flipud(circshift([X],3,2)), 'PolarDirection', 'ccw','AngularRange', [0 2*pi]-pi/12,'TickSpacing', 360,'GridStyle', '-', 'PolarGrid', {0 0},'MeshScale', [1 1]-1./size(X));
if isfield(cfg, 'clim')
  CL = cfg.clim;
  clim(CL);
elseif sign(min(X(:))) ~= sign(max(X(:)))
  CL = [-1 1]* max(abs(X(:)));
  clim(CL);
else
  CL = clim;
end
view([0,90])
axis off
box off
% colormap((brewermap(256, '')))
colormap(inferno)
%% colorbar
ax(4) = axes('Position', [.65 .1 .3 .8]);
tmp=imagesc(X, CL);
cb = colorbar;
tmp.Visible='off';
box off, axis off
cb.Label.String = cfg.cblabel;
cb.Label.FontSize = 12;
cb.FontSize = 10;
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
ax(3).FontSize=10;
% ax(3).Color = 'w';
set(gca,'XColor',[1 1 1]);
xlabel(cfg.timelabel)


