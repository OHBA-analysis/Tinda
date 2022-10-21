function ax = cyclicalstateplot_perstate(ordering,mean_direction,sigpoints,plotstates,newfigure,color_scheme)
% Plot state network as circular diagram with arrows

if nargin<1
  error('Need ordering of states to plot')
end
K = length(ordering);
if ~exist('plotstates', 'var') || isempty(plotstates)
  plotstates=1:K;
  row=3; col=4;
end
if exist('newfigure','var') && newfigure==true || ~exist('newfigure','var') 
  if length(plotstates)>1
    setup_figure([],2,1);
  else
    figure('Position',[440 501 402 297]);
  end
else
 axes(gca);
end

if nargin<6
    color_scheme = set1_cols;
end
disttoplot_manual = zeros(12,2);
for i=1:12
  temp = exp(sqrt(-1)*(i+2)/12*2*pi);
  disttoplot_manual(ordering(i),:) = [real(temp),imag(temp)];
end

% ordering = [ordering(1),fliplr(ordering(2:end))]; % ensures rotates in clockwise direction

for k_to_plot = plotstates
  if length(plotstates)>1
    subplot(row,col,k_to_plot);
  end
  for ik1=ordering(k_to_plot)
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
            linescale*[disttoplot_manual(k2,1)-disttoplot_manual(ik1,1)],linescale*[disttoplot_manual(k2,2)-disttoplot_manual(ik1,2)],'Color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.4/quivlength);hold on;
        else % arrow from k2 to k1:
          quiver(disttoplot_manual(k2,1),disttoplot_manual(k2,2),...
            linescale*[disttoplot_manual(ik1,1)-disttoplot_manual(k2,1)],linescale*[disttoplot_manual(ik1,2)-disttoplot_manual(k2,2)],'Color',[0 0 0]+0.8,'LineWidth',2,'MaxHeadSize',0.4/quivlength);hold on;
        end
      end
    end
  end
  
  % do it once more for state k in the columns (switch around the arrow
  % colors)
  for ik1=1:K
    for k2=ordering(k_to_plot)
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
            linescale*[disttoplot_manual(k2,1)-disttoplot_manual(ik1,1)],linescale*[disttoplot_manual(k2,2)-disttoplot_manual(ik1,2)],'Color',[0 0 0]+0.8,'LineWidth',2,'MaxHeadSize',0.4/quivlength);hold on;
        else % arrow from k2 to k1:
          quiver(disttoplot_manual(k2,1),disttoplot_manual(k2,2),...
            linescale*[disttoplot_manual(ik1,1)-disttoplot_manual(k2,1)],linescale*[disttoplot_manual(ik1,2)-disttoplot_manual(k2,2)],'Color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.4/quivlength);hold on;
        end
      end
    end
  end

  % plot the circles with numbers
  for ik=1:K
    docolor=1;
    if docolor, edgecol = color_scheme{ik}; facecol = edgecol; else, edgecol = [0 0 0]; facecol = [1 1 1]; end
    scatter1(ik) = scatter(disttoplot_manual(ik,1),disttoplot_manual(ik,2),400,...
      'MarkerFaceColor', facecol,'MarkerEdgeColor', edgecol);
    hold on
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatter1(ik).MarkerFaceAlpha = 1;%.75;
    t = text(disttoplot_manual(ik,1),disttoplot_manual(ik,2),int2str((ik)),'FontSize',12,'FontWeight','bold', 'HorizontalAlignment', 'center', 'Visible', 'On', 'FontName', 'Calibri');
end
  axis square
  axis off
end
end