
whichstudy=6;
load(sprintf('/Users/matsvanes/Documents/Werk/Tinda/study%d.mat', whichstudy))
addpath(genpath('/Users/matsvanes/Werk/scripts/Tinda/'))
color_scheme = colorscheme(whichstudy);
studies = {'MEG UK', '', 'HCP', '', '', 'Cam-CAN'};

if whichstudy==1
  numBins=8;
elseif whichstudy==3
  numBins=8;
elseif whichstudy==6
  numBins=16;
end

do_metastate=1;
%%
for do_image=[true, false]
  for do_demean_FO=true
    FO=hmm_1stlevel.FO;
    if do_demean_FO
      FO = FO-mean(FO);
      ttl1 = {"Subjects' circle density relative to group mean", ""};
    else
      ttl1 = {"Subjects' circle density", ''};
    end
    
    K=12;
    theta = (0:1/12:.99)*2*pi; % define angles
    theta = circshift(theta,-3); % set first angle to 12 o'clock
    for k=1:K % reshuffle angles in bestseq order
      tmp(k) = theta(bestseq==k);
    end
    theta = tmp; % theta now corresponds to the angle of each state in the cycle plot
    x = cos(theta); % its x-value
    y = sin(theta); % its y-vale

    % find the density for each subject
    x_weighted = FO.*x;
    y_weighted = FO.*y;

    x_weighted_average = mean(x_weighted,2);
    y_weighted_average = mean(y_weighted,2);

    %% circle density all subjects

    [xCenters, yCenters, counts] = histcounts2_circle(x_weighted_average, y_weighted_average, numBins);

    fig=setup_figure([], 2,1);
    ax(1) = axes('Position', [.1 .1 .8 .8]);
    if do_image % Display the density map using imagesc
      xpos = .06;
      imagesc(xCenters, yCenters, counts');
      axis xy;   % To ensure y-axis is in the correct direction
      colorbar;  % To show color scale
    else
      xpos=.1;
      polarhistogram(angle(x_weighted_average+sqrt(-1)*y_weighted_average), -pi/12:pi/6:11*pi/6);
      thetaticks([])
    end
    ax(2) = axes('Position', [xpos 0.1 .8 .8]);
    cyclicalstateplot(bestseq, hmm_1stlevel.cycle_metrics.mean_assym, zeros(12), color_scheme)
    title(ttl1)


    fname = sprintf('/Users/matsvanes/Documents/Werk/Tinda/figures/circle_density_study%d', whichstudy);
    if do_demean_FO
      fname = [fname, '_demean'];
    end
    if do_image
      fname = [fname, '_image'];
    else
      fname = [fname, '_rose'];
    end
%     print(fname, '-dpng')

    %% Circle density difference between groups
    %
    %strongest vs weakest rotational momentum
    [val, group1] = sort(hmm_1stlevel.cycle_metrics.rotational_momentum, 'ascend');
    [val, group2] = sort(hmm_1stlevel.cycle_metrics.rotational_momentum, 'descend');

    if whichstudy==1
      group1=group1(1:27);
      group2=group2(1:27);
      if do_demean_FO
        ttl2 = {'Circle density relative to mean:', '27 strongest-weakest M', ''};
      else
        ttl2={'Circle density: 27 strongest minus 27 weakest M', ''};
      end
      lg = {'Weakest M', 'Strongest M'};
    elseif whichstudy==3

      if do_metastate
        group1=clustermember==0; %cognitive
        group2=clustermember==1; %sensorimotor
        if do_demean_FO
          ttl2 = {'Circle density relative to mean:', 'cognitive-sensorimotor metastate subjects', ''};
        else
          ttl2 = {'Circle density: cognitive-sensorimotor metastate subjects', ''};
        end
        lg = {'Sensorimotor', 'Cognitive'};
      else
        group1=group1(1:39);
        group2=group2(1:39);
        if do_demean_FO
          ttl2 = {'Circle density relative to mean:', '39 strongest-weakest M', ''};
        else
          ttl2={'Circle density: 39 strongest minus 39 weakest M', ''};
        end
        lg = {'Weakest M', 'Strongest M'};
      end
    elseif whichstudy==6
      group1=group1(1:150);
      group2=group2(1:150);
      group1 = info(:,1)<40;
      group2 = info(:,1)>60;
      if do_demean_FO
        ttl2 = {'Circle density relative to mean:', '150 strongest-weakest M', ''};
      else
        ttl2={'Circle density: 150 strongest minus 150 weakest M', ''};
      end
      lg = {'Weakest M', 'Strongest M'};
    end
        % test for nonuniformity

        diff = [x_weighted_average(group1)+sqrt(-1)*y_weighted_average(group1); -(x_weighted_average(group2)+sqrt(-1)*y_weighted_average(group2))];
        [pval, z] = circ_rtest(angle(diff));


    [xCenters, yCenters, counts1] = histcounts2_circle(x_weighted_average(group1), y_weighted_average(group1), numBins);
    [xCenters, yCenters, counts2] = histcounts2_circle(x_weighted_average(group2), y_weighted_average(group2), numBins);
    counts = counts1-counts2;

    crange = .8*max(abs(counts(:)))*[-1 1];

    fig=setup_figure([], 2,1);
    ax(1) = axes('Position', [.1 .1 .8 .8]);
    if do_image  % Display the density map using imagesc
      xpos=.06;
      imagesc(xCenters, yCenters, counts', crange);
      colormap(flipud(brewermap(64,'RdBu')))
      axis xy;   % To ensure y-axis is in the correct direction
      colorbar;  % To show color scale
    else
      xpos=.1;
      polarhistogram(angle(x_weighted_average(group2)+sqrt(-1)*y_weighted_average(group2)),-pi/12:pi/6:11*pi/6);
      hold on
      polarhistogram(angle(x_weighted_average(group1)+sqrt(-1)*y_weighted_average(group1)),-pi/12:pi/6:11*pi/6);
      thetaticks([])
      ax(3)=axes('Position', [.95,.95, .001,.001]), plot([1:2; 1:2],[1:2;1:2]), axis off, box off, xticks([]), yticks([]), legend(lg)
    end

    ax(2) = axes('Position', [xpos 0.1 .8 .8]);
    cyclicalstateplot(bestseq, hmm_1stlevel.cycle_metrics.mean_assym, zeros(12), color_scheme)
    title(ttl2)

    if do_metastate
      fname = sprintf('/Users/matsvanes/Documents/Werk/Tinda/figures/circle_density_diff_metastate_study%d', whichstudy);
    else
      fname = sprintf('/Users/matsvanes/Documents/Werk/Tinda/figures/circle_density_diff_study%d', whichstudy);
    end
    if do_demean_FO
      fname = [fname, '_demean'];
    end
    if do_image
      fname = [fname, '_image'];
    else
      fname = [fname, '_rose'];
    end
%     print(fname, '-dpng')

  end
end




%%
FO=hmm_1stlevel.FO;
[h,p,ci,stats] = ttest2(FO(clustermember==1,:), FO(clustermember==0,:), 'Alpha',0.05/12);

p_thresh = p;
p_thresh(h==0) = 1;


fig=setup_figure([], 2,1);
y=circshift(sin((0:1/12:.99)*2*pi)/2.5+0.425, -3);
x=circshift(cos((0:1/12:.99)*2*pi)/2.5+0.425, -3);
for k=1:12
% subplot(3,4,k), 
% boxplot_with_scatter(FO(:,k), [],[],2-clustermember); title(sprintf('state %d',k)); 
% sigstar([1 2], p_thresh(k));
ax(k) = axes('Position', [x(k), y(k) 0.15 0.15]);
boxplot_with_scatter(FO(:,bestseq(k)), [],[],2-clustermember); title(sprintf('state %d',bestseq(k))); 
sigstar([1 2], p_thresh(bestseq(k)));

xticklabels({'Sens', 'Cog'})
end
%%
suptitle('FO for sensori/motor (blue) and cognitive (red) metastate subjects')
if do_metastate && whichstudy==3
  save_figure(fig, sprintf('/Users/matsvanes/Documents/Werk/Tinda/figures/FO_diff_metastate_study%d', whichstudy), false)
else
  save_figure(fig, sprintf('/Users/matsvanes/Documents/Werk/Tinda/figures/FO_diff_study%d', whichstudy), false)
end
%%
% clustering of the FO correlation matrix
FO = hmm_1stlevel.FO;
FOcorr = corr(FO);
[A, B] = nets_hierarchy(FOcorr, FOcorr, 1:12, '');

% b(:,1) is now each subject's strength of clustering
[a,b] = pca(FO,'NumComponents',2);
clustermember_meg = b(:,1)>0;


gmmfit = fitgmdist(b,2);
prob_fit = posterior(gmmfit,b);
prob_fit = prob_fit(:,2);
clustermember_meg = prob_fit<=0.5; 