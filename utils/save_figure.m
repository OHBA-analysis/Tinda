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

% create directory if it doesn't exist
if ~exist(fileparts(path), 'dir')
  mkdir(fileparts(path))
end


print(path,'-depsc')
print(path,'-dpng')
print(path,'-dsvg')
if highres
  set(fig, 'Position', [1 1 2 2].*fig.Position)
  set_font(14);
  print([path '_highres'],'-depsc');
  print([path '_highres'],'-dpng');
end