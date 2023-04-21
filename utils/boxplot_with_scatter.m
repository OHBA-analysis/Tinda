function boxplot_with_scatter(data, facecolor, facealpha, group)

if ~exist('facealpha', 'var') || isempty(facealpha)
    facealpha=0.6;
end

if ~exist('group', 'var'),
    group=[];
end

[n, m] = size(data);
hold on
if m==1
    if ~isempty(group)
        ug = unique(group)
        for k2=1:length(ug)
            k3=ug(k2)
            if ~exist('facecolor', 'var') || isempty(facecolor)
                scatter(1+ones(sum(group==k3),1).*((k3-1)+(rand(sum(group==k3),1)-0.5)/2),data(group==k3,1),'filled', 'MarkerFaceAlpha',facealpha)
            else
                if iscell(facecolor)
                    MarkerFaceColor = facecolor{k3};
                else
                    MarkerFaceColor = facecolor;
                end
                scatter(1+ones(sum(group==k3),1).*((k3-1)+(rand(sum(group==k3),1)-0.5)/2),data(group==k3,1),'filled','MarkerFaceColor', MarkerFaceColor, 'MarkerFaceAlpha',facealpha)
            end
        end
    else
        k=1;scatter(1+ones(size(data(:,k))).*((k-1)+(rand(size(data(:,k)))-0.5)/2),data(:,k),'filled', 'MarkerFaceAlpha',facealpha)
    end
else
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
end
if exist('group', 'var') && ~isempty(group)
    h = boxplot(data, group, 'Colors', 'k', 'Width', .75, 'whisker', 1000); % put this on top
else
    h = boxplot(data, 'Colors', 'k', 'Width', .75, 'whisker', 1000); % put this on top
end
set(h, {'linew'}, {2})