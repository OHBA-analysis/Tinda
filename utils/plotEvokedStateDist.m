function fig1 = plotEvokedStateDist(gam_evoked_sj,t)

K = size(gam_evoked_sj,2);
color_scheme = set1_cols;
meanresponse = mean(gam_evoked_sj,3);
ste = std(gam_evoked_sj,[],3) ./ sqrt(size(gam_evoked_sj,3));
figure();

ls = {'-','--'};
for k=1:K
    shadedErrorBar(t,meanresponse(:,k),ste(:,k),{'LineWidth',2','Color',color_scheme{k},'LineStyle',ls{1+(k>6)}},0.5);
    hold on;
    h(k) = plot(NaN,NaN,'Color',color_scheme{k},'LineWidth',2,'LineStyle',ls{1+(k>6)});
end
%xlim([1,2.5]);
ylim([-0.15,0.15]);
% do cluster sign flipping test:
clear corrp tstats_cluster
for k=1:K
    dat = permute(gam_evoked_sj(:,k,:),[3,2,1]);
    thresh = 3;
    nP=5000;
    [corrp(:,k),tstats_cluster(:,k)] = osl_clustertf(dat,thresh,nP);
end
hyp = corrp>(1-1e-3);
ylow = ylim;

ylow(1) = -0.0689;

yvals = ylow(1) * [1.05:0.05:2];
lengththresh = 2; % do not plot clusters less than this length
for iSt=1:K
    sig_line = find(diff([0;hyp(:,iSt);0]));
    S = reshape(sig_line,[2,length(sig_line)/2]);
    S(:,diff(S)<lengththresh) = [];
    sig_line = S(:);
    sig_line(sig_line>length(hyp))=length(hyp);
    yval = yvals(iSt);
    for k=1:length(sig_line)/2
        line([t(sig_line(k*2-1)),t(sig_line(k*2))],[yval,yval],'Color',color_scheme{iSt},'LineWidth',2);
    end
end
ylim([yval*1.05,ylow(2)]);
for k=1:K,h(k).DisplayName=['RSN-State ',int2str(k)];end
leg=legend(h,'Location','EastOutside');
plot4paper('Time (sec)','Relative state probability');

end