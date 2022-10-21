function fig = cyclicalstate_distributionplot(ordering, weights, CL)

% edges = logspace(log10(CL(1)), log10(CL(2)), 100);
edges = linspace(CL(1), 1.1*CL(2), 100);
color_ix = discretize(weights, edges);

CM = colormap(inferno(120));
CM = CM(20:end,:);
for k=1:12
  color_scheme{k} = CM(color_ix(k),:);
end
cyclicalstatesubplot(ordering, zeros(12), zeros(12), color_scheme)
hold on
spider_plot(color_ix(ordering)', 'AxesHandle', gca,'AxesLimits', repmat([1; 150], [1,12]), 'AxesLabels', 'none', 'Direction', 'counterclockwise', 'FillOption', 'on', 'AxesDisplay','none', 'Marker', 'none', 'AxesColor',[1 1 1])
