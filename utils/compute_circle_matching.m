studies = [1,3,6,7];
studylabel{1} = 'MEG UK';
studylabel{3} = 'HCP';
studylabel{6} = 'Cam-CAN';
studylabel{7} = 'Wakeman-Henson';
studylabel{8} = 'Wakeman-Henson';

template_study=1;
for istudy=studies
    cf{istudy} = getStudyDetails(istudy);
    map{istudy} = load([cf{istudy}.resultsdir, 'group_avg_spectral_maps.mat']);
    if istudy==1
      bs{istudy} = load([cf{istudy}.resultsdir, sprintf('bestseq%d_coherence.mat', istudy)]);
    else
      bs{istudy} = load([cf{istudy}.resultsdir, sprintf('bestseq%d_study1matched', istudy)]);
    end
    bs{istudy} = bs{istudy}.bestsequencemetrics{1};
    tmp = circle_angles(bs{istudy});
    angles{istudy} = tmp(:,1);
end


nperm=10000;
for istudy=setdiff(studies, template_study)
    delta = angles{istudy}-angles{template_study};

    obs = mean(abs(delta));
    stat_bestseq_similarity{istudy}.obs = obs;
    
    perm=zeros(nperm,1);
    for iperm=1:nperm
        tmp = circle_angles(randperm(12));
        perm(iperm,1) = mean(abs(angles{template_study}-tmp(:,1)));
    end
    stat_bestseq_similarity{istudy}.perm = perm;
    stat_bestseq_similarity{istudy}.pval = max([sum(perm<obs)/nperm, 1/(nperm+1)]);
end