function save_figure(fig, path, highres)
% this function saves a figure in two formats, in two sizes.
if ~exist('fig', 'var') || isempty(fig)
  fig = gcf;
elseif ischar(fig)
  if exist('path', 'var')
    highres=path;
  end
  path = fig;
  fig = gcf;
end
if ~exist('highres', 'var'), highres=true; end

% 
print(path,'-depsc')
print(path,'-dpng')
if highres
  set(fig, 'Position', [1 1 2 2].*fig.Position)
  print([path '_highres'],'-depsc');
  print([path '_highres'],'-dpng');
end