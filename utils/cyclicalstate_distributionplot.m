function fig = cyclicalstate_distributionplot(ordering, weights, CL, CM, cycleplot, spider_color, lim, offset)

% edges = logspace(log10(CL(1)), log10(CL(2)), 100);
edges = linspace(CL(1), 1.1*CL(2), 100);
color_ix = discretize(weights, edges);

if ~exist('CM', 'var') || isempty(CM)
  CM = colormap(inferno(120));
  CM = CM(20:end,:);
end
for k=1:12
  color_scheme{k} = CM(color_ix(k),:);
end
if ~exist('cycleplot') || isempty(cycleplot) || cycleplot==true
  cyclicalstatesubplot(ordering, zeros(12), zeros(12), color_scheme)
end
hold on

if ~exist('spider_color') || isempty(spider_color)
  clr = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};
  spider_color = clr{1};
end
if ~exist('lim', 'var') || isempty(lim)
    lim = length(CM);
end
if ~exist('offset', 'var') || isempty(offset)
    offset = 1;
end
spider_plot(color_ix(ordering)', 'AxesHandle', gca, 'AxesLimits', repmat([1; lim], [1,12]), 'Color', spider_color,'AxesLabels','none', 'AxesOffset', offset, 'Direction', 'counterclockwise', 'FillOption', 'on', 'AxesDisplay','none', 'Marker', 'none', 'AxesColor',[1 1 1])

