function ax = cyclicalstateplot(ordering,mean_direction,sigpoints,color_scheme,newfigure,plotstate,allowstretch, arrowcol)
% Plot state network as circular diagram with arrows

if nargin<1
    error('Need ordering of states to plot')
end
K = length(ordering);
clock=[12,1,2,3,11,4,10,5,7,9,6,8];
if nargin<4 || isempty(color_scheme)
    color_scheme = set1_cols;
end
if nargin<5,  newfigure=true;   end
if ~exist('arrowcol', 'var') || isempty(arrowcol)
  arrowcol = [0 0 0] ; % or [0.8 0.8 0.8] for grey
end
arrowcol2 = arrowcol;
if exist('plotstate', 'var') && length(plotstate)==1
  arrowcol2=arrowcol+0.8;
else
  plotstate=1:K;
end
if ~exist('allowstretch', 'var'), allowstretch=0; end

disttoplot_manual = zeros(K,2);
for i=1:K
    temp = exp(sqrt(-1)*(i+2)/K*2*pi);
    disttoplot_manual(ordering(i),:) = [real(temp),imag(temp)];
end

if newfigure
 figure('Position',[440 501 402 297]);
else
 axes(gca);
end
tmpa = gca;


for ik1=plotstate
    for k2=1:K
        if sigpoints(ik1,k2)
            linescale = sqrt(sum((disttoplot_manual(k2,:)-disttoplot_manual(ik1,:)).^2));
            if linescale>1.42 % ie line is four or more steps
                linescale = 1;
            elseif linescale >1.2% line is three steps
                linescale = 0.98;
            elseif linescale > 0.75 % line is two steps
                linescale = 0.9;
            else  % line is one step
                linescale = 0.8;
            end
            quivlength = linescale * sqrt(sum((disttoplot_manual(k2,:)-disttoplot_manual(ik1,:)).^2));
            if mean_direction(ik1,k2)>0 % arrow from k1 to k2:
                  quiver([disttoplot_manual(ik1,1)],[disttoplot_manual(ik1,2),],...
                      linescale*[disttoplot_manual(k2,1)-disttoplot_manual(ik1,1)],linescale*[disttoplot_manual(k2,2)-disttoplot_manual(ik1,2)],'Color',arrowcol,'LineWidth',2,'MaxHeadSize',0.4/quivlength);hold on;
            else % arrow from k2 to k1:
                quiver(disttoplot_manual(k2,1),disttoplot_manual(k2,2),...
                    linescale*[disttoplot_manual(ik1,1)-disttoplot_manual(k2,1)],linescale*[disttoplot_manual(ik1,2)-disttoplot_manual(k2,2)],'Color',arrowcol2,'LineWidth',2,'MaxHeadSize',0.4/quivlength);hold on;
            end
        end
    end
end

msize = (0.5+tmpa.Position(3)/2)*400;
for ik=1:K
  docolor=1;
  if docolor, edgecol = color_scheme{ik}; facecol = edgecol; else, edgecol = [0 0 0]; facecol = [1 1 1]; end
    scatter1 = scatter(disttoplot_manual(ik,1),disttoplot_manual(ik,2),msize,...
        'MarkerFaceColor', facecol,'MarkerEdgeColor', edgecol); 
    hold on
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatter1.MarkerFaceAlpha = 1;%.75;
    text(disttoplot_manual(ik,1),disttoplot_manual(ik,2),int2str((ik)),'FontSize',12,'FontWeight','bold', 'FontName', 'Calibri', 'HorizontalAlignment', 'center');hold on;
end
if ~allowstretch
  axis square
end
axis off
end