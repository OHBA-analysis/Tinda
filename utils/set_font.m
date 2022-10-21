function fig = set_font(fig, sz, foc)
if nargin<1 || isempty(fig)
  fig = gcf;
elseif nargin<2 && isnumeric(fig)
  sz = fig;
  fig = gcf;
elseif nargin<3 && isnumeric(fig)
  foc = sz;
  sz = fig;
  fig = gcf;
end
if ~exist('sz','var') || isempty(sz)
  sz = 10;
end
if ~exist('foc', 'var') || isempty(foc)
  foc = 'all';
end

if ischar(foc)
  foc = {foc};
end
for kk=1:length(foc)
  f = foc{kk};
  switch f
    case 'all'
      set(findall(fig,'-property','FontSize'),'FontSize',sz) % adjust fontsize to your document
    case 'label'
      for k=1:length(fig.Children)
        try
          fig.Children(k).XLabel.FontSize = sz;
        end
        try
          fig.Children(k).YLabel.FontSize = sz;
        end
      end
    case 'title'
      for k=1:length(fig.Children)
        try
          fig.Children(k).Title.FontSize = sz;
        end
      end
  end
end
set(findall(fig,'-property','FontName'),'FontName','Calibri')
