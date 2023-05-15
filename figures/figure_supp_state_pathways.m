
use_WB_nnmf=true; % whether or not to use the wide band NNMF Diego's Nature Comms to select power and coherence (alternative is selecting 1-30 Hz)
useMT = true; %

fig=setup_figure([],2,1.5);
clock = [12, 1:11];
cnt1=1;
studies = {'MEG UK', 'Cam-CAN', 'HCP'};
for whichstudy=[1,6,3]
    config = getStudyDetails(whichstudy);
    color_scheme = colorscheme(whichstudy);
    load(config.metricfile, 'hmm_1stlevel')
    
    if strcmp(config.reordering_states, 'coherence')
        optimalseqfile = [config.resultsdir,'bestseq',int2str(whichstudy),'_coherence' ,'.mat'];
    else
        optimalseqfile = [config.hmmfolder,'bestseq',int2str(whichstudy),'.mat'];
    end
    load(optimalseqfile);
    bestseq = bestsequencemetrics{1};
    angleplot = circle_angles(bestseq);
    
    bestseq_clockwise = [bestseq(1), bestseq(end:-1:2)];
       
    K=length(bestseq);
    sigpoints =  hmm_1stlevel.assym_ttest.sigpoints;

    for k=1:K
        if k<=6
            ax(k, cnt1) = axes('Position', [0.01+(cnt1-1)/6, .85-.9*((k-1)/(K/2)), 0.85*1/6, 1.1*1/K]);
        else
            ax(k, cnt1) = axes('Position', [0.51+(cnt1-1)/6, .85-.9*(((k-1)-K/2)/(K/2)), 0.85*1/6, 1.1*1/K]);
        end
        mask = zeros(K);
        mask(bestseq_clockwise(k), :) = 1;
        mask(:, bestseq_clockwise(k)) = 1;
        cyclicalstateplot(bestseq, hmm_1stlevel.cycle_metrics.mean_direction, sigpoints.*mask, color_scheme, false,[],false,[],200);
        if whichstudy==6
            title({sprintf("%d o'clock", clock(k)), ''})
        end
    end
    
    
    
    ax(K+1, cnt1) = axes('Position', [0.01+(cnt1-1)/6, .89, 0.85*1/6, 1.1*1/K]);
    axis off
    title(sprintf('%s', studies{cnt1}))
    ax(K+2, cnt1) = axes('Position', [0.51+(cnt1-1)/6, .89, 0.85*1/6, 1.1*1/K]);
    axis off
    title(sprintf('%s', studies{cnt1}))
    
    cnt1=cnt1+1;
end
for whichstudy=[1,6,3]
    config = getStudyDetails(whichstudy);
    save_figure([config.figdir, 'figure1_tinda_example/', '1supp_StatePathways_allstudies'],[],false);
end