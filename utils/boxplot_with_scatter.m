function boxplot_with_scatter(data, facecolor, facealpha)

if ~exist('facealpha', 'var') || isempty(facealpha)
  facealpha=0.6;
end

[n, m] = size(data);
hold on
for k=1:m
  if ~exist('facecolor', 'var') || isempty(facecolor)
    scatter(1+ones(size(data(:,k))).*((k-1)+(rand(size(data(:,k)))-0.5)/2),data(:,k),'filled', 'MarkerFaceAlpha',facealpha)
  else
    if iscell(facecolor)
      MarkerFaceColor = facecolor{k};
    else
      MarkerFaceColor = facecolor;
    end
    scatter(1+ones(size(data(:,k))).*((k-1)+(rand(size(data(:,k)))-0.5)/2),data(:,k),'filled', 'MarkerFaceAlpha',facealpha, 'MarkerFaceColor',MarkerFaceColor)
  end
end
h = boxplot(data, 'Colors', 'k', 'Width', .75, 'whisker', 1000); % put this on top
set(h, {'linew'}, {2})