function ax = cyclicalstateplot(ordering,mean_direction,sigpoints,color_scheme)
% Plot state network as circular diagram with arrows

if nargin<1
    error('Need ordering of states to plot')
end
K = length(ordering);
clock=[12,1,2,3,11,4,10,5,7,9,6,8];
if nargin<4
    color_scheme = set1_cols;
end
disttoplot_manual = zeros(12,2);
for i=1:12
    temp = exp(sqrt(-1)*(i+2)/12*2*pi);
    disttoplot_manual(ordering(i),:) = [real(temp),imag(temp)];
end


figure('Position',[440 501 402 297]);
for ik1=1:K
    for k2=1:K
        if sigpoints(ik1,k2)
%             line([disttoplot(ik,1),disttoplot(jk,1)],[disttoplot(ik,2),disttoplot(jk,2)],...
%                 'color',0.5*[1,1,1]);hold on;
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
                    linescale*[disttoplot_manual(ik1,1)-disttoplot_manual(k2,1)],linescale*[disttoplot_manual(ik1,2)-disttoplot_manual(k2,2)],'Color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.4/quivlength);hold on;
            end
        end
    end
end
for ik=1:K
%     scatter1 = scatter(disttoplot_manual(ik,1),disttoplot_manual(ik,2),400,...
%         'MarkerFaceColor',color_scheme{ik},'MarkerEdgeColor',color_scheme{ik}); 
    scatter1 = scatter(disttoplot_manual(ik,1),disttoplot_manual(ik,2),400,...
        'MarkerFaceColor',[1 1 1],'MarkerEdgeColor','k'); 
    hold on
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatter1.MarkerFaceAlpha = 1;%.75;
    if ik==10
%         scatter1 = scatter(disttoplot_manual(ik-1,1),disttoplot_manual(ik-1,2),400,...
%         'MarkerFaceColor',color_scheme{ik-1},'MarkerEdgeColor',color_scheme{ik-1}); 
    scatter1 = scatter(disttoplot_manual(ik,1),disttoplot_manual(ik,2),400,...
        'MarkerFaceColor',[1 1 1],'MarkerEdgeColor','k');         
    hold on
        % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
        scatter1.MarkerFaceAlpha = 0.5;
        text(disttoplot_manual(ik-1,1)-0.03,disttoplot_manual(ik-1,2),int2str(ik-1),'FontSize',12,'FontWeight','bold');hold on;
    end
    if ik<10
%         text(disttoplot_manual(ik,1)-0.03,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
        text(disttoplot_manual(ik,1)-0.03,disttoplot_manual(ik,2),int2str(clock(ik)),'FontSize',12,'FontWeight','bold');hold on;
    else
%         text(disttoplot_manual(ik,1)-0.05,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
text(disttoplot_manual(ik,1)-0.05,disttoplot_manual(ik,2),int2str(clock(ik)),'FontSize',12,'FontWeight','bold');hold on;
    end
end
axis square
axis off
end