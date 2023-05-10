% load('/Users/matsvanes/Downloads/replay_bubblechart.mat')
load('/Users/matsvanes/Downloads/tinda_replay_perc5_2.mat')
addpath('utils/')
addpath('/Users/matsvanes/Documents/Werk/scripts/sigstar')
whichstudy=1;
color_scheme = colorscheme(whichstudy);


for k=1:2
  assym{k} = transpose(squeeze(replay.K13.study{k}.FO(13,1:12,1,:)-replay.K13.study{k}.FO(13,1:12,2,:)));
end

%% these are the percentiled versions



for d=1

  for istudy=1:2
    fig = setup_figure([],2,.8);
    
    for ints=1:5
      IT(ints) = round(mean(squash(cellfun(@mean, replay.K13.perc.t_intervals{ints}))));
    end


    ax(istudy) = axes('Position', [0.075 .1, .2, .8]);
    boxplot_with_scatter(assym{istudy}, color_scheme);
    for k=1:12
      group{k} = [k k];
    end
    sigstar(group, replay.K13.study{istudy}.assym_ttest.pvals(13,1:12)*24);

    hline(0)
    xlim([.1 12.9])
    box off
    view([90,90])
    xlabel('State')
    %   title({'State FO assymmetry', 'in replay intervals'})
    ylabel('FO asymmetry')
    if istudy==1
      xlabel({'\bf \fontsize{15} Study 1 \rm', '\fontsize{12} State'})
      ylim([-0.025 .03])
    else
      xlabel({'\bf \fontsize{15} Study 2 \rm', '\fontsize{12} State'})
      ylim([-0.025 .045])
    end
    text(0.9, 1.075, 'State FO asymmetry in replay intervals', 'Units', 'normalized', 'FontSize', 14, 'FontWeight','bold')

    %
    ax(istudy+3) = axes('Position', [.25 .15 .8 .75]);
    cyclicalstateplot(replay.K13.perc.bubbleplot.bestseq_nott, zeros(12), zeros(12), color_scheme, false)
    ax(istudy+2) = axes('Position', [.25 .175 .8 .7]);

    cmap = replay.K13.perc.bubbleplot.cmap;
    % deal with nans by appending a white to the end
    cmap = [cmap; 1,1,1];
    bestseq = replay.K13.perc.bubbleplot.bestseq_nott;

    %   metric = replay.K13.perc.bubbleplot.tstat(:,:,2:end);
    metric = replay.K13.perc.study{istudy}.mean_assym;
    % metric = replay.K13.perc.bubbleplot.mean_direction;
    th = angle(replay.K13.perc.bubbleplot.angleplot_nott(:,1));
    th=(2*pi - th) + 0.25*2*pi; % reconfigure



    if d==1 % positive value means replay comes after
      metric_sel = squeeze(metric(13,1:12,:));
    else
      metric_sel = squeeze(metric(1:12,13,:));
    end
    ls = linspace(-max(abs(squash(metric_sel))),max(abs(squash(metric_sel))),256);
    metric_discrete=discretize(metric_sel,ls);

    metric_discrete(isnan(metric_discrete)) = 257;


    pr=10;
    ir=10;
    isz = .7;
    psz = 1;
    sz_exp = .8;


    k=0;
    C0 = repmat([1,1,1], [12,1]);
    R0 = ir+pr*k*ones(12,1);
    SZ0 = ones(12,1);

    k=1;
    C1 = cmap(metric_discrete(:,k),:);
    R1 = ir+pr*k*ones(12,1);
    SZ1 = isz+psz*k^sz_exp*ones(12,1);

    k=2;
    C2 = cmap(metric_discrete(:,k),:);
    R2 = ir+pr*k*ones(12,1);
    SZ2 = isz+psz*k^sz_exp*ones(12,1);

    k=3;
    C3 = cmap(metric_discrete(:,k),:);
    R3 = ir+pr*k*ones(12,1);
    SZ3 = isz+psz*k^sz_exp*ones(12,1);

    k=4;
    C4 = cmap(metric_discrete(:,k),:);
    R4 = ir+pr*k*ones(12,1);
    SZ4 = isz+psz*k^sz_exp*ones(12,1);

    k=5;
    C5 = cmap(metric_discrete(:,k),:);
    R5 = ir+pr*k*ones(12,1);
    SZ5 = isz+psz*k^sz_exp*ones(12,1);

    tbl = table(th, C0,R0,SZ0, C1,R1,SZ1, C2,R2,SZ2, C3,R3,SZ3, C4,R4,SZ4, C5,R5,SZ5);

    if d==2
      fig = setup_figure([],2,1);
      ax(1) = axes('Position', [.05 .0 .8 .9]);
    end

    h=polarbubblechart(tbl,'th',{'R0', 'R1', 'R2', 'R3', 'R4', 'R5'},{'SZ0', 'SZ1', 'SZ2','SZ3','SZ4', 'SZ5'}, {'C0', 'C1', 'C2', 'C3', 'C4', 'C5'}, 'MarkerFaceAlpha',1);% thetaticklabels(circshift(bestseq,3))
    thetaticklabels([])
    fig.CurrentAxes.ThetaAxis.FontWeight='bold';
    fig.CurrentAxes.ThetaAxis.FontSize = 14;
    rlim([0 60])
    rticks(linspace(20,60,5))
    rticklabels(IT)
    fig.CurrentAxes.RAxisLocation=0;
    ax(2) = axes('Position', [.4 .0 .5 0.075]);
    colormap(cmap);
    ax(2).CLim = [-1,1]* max(abs(squash(metric_sel)));
    cb = colorbar('Location','North');
    % cb.Position(3) = cb.Position(3)*3;
    cb.Position(4) = cb.Position(4)*2;
    cb.Title.String = 'FO asymmetry';
    cb.Title.FontSize=14;
    cb.Title.FontWeight = 'bold';
    axis off
    cb.Label.FontSize=14;
    cb.Label.FontWeight = 'bold';
    if d==1
      fname = '/Users/matsvanes/Documents/Werk/Tinda/Study1/figures/replay/bubbleplot_replay_interval_study';
    else
      fname = '/Users/matsvanes/Documents/Werk/Tinda/Study1/figures/replay/bubbleplot_state_interval';
    end

    save_figure([fname num2str(istudy)],[],false)
  end
end