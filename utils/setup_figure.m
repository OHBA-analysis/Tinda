function fig = setup_figure(fig,width,hw_ratio,box)
if ~exist('fig', 'var') || isempty(fig)
  fig=figure;
end
if ~exist('width', 'var') || isempty(width)
  width=2; % columns
end
if ~exist('hw_ratio', 'var') || isempty(hw_ratio)
  hw_ratio=0.5; 
end
if ~exist('box', 'var') || isempty(box)
  box='off'; 
end

% Elsevier recommendations
if width==1
  picturewidth = 9; 
elseif width==1.5
  picturewidth = 14; 
elseif width==2
  picturewidth = 19; 
end
set_font(10)

set(findall(fig,'-property','Box'),'Box',box) % optional
set(findall(fig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(fig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
