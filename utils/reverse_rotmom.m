for whichstudy=[3 4 6]
    keep whichstudy
    config = getStudyDetails(whichstudy);
    load(config.metricfile, 'hmm_1stlevel')
    % basic
    hmm_1stlevel.cycle_metrics.rotational_momentum = -hmm_1stlevel.cycle_metrics.rotational_momentum;
    hmm_1stlevel.cycle_metrics.rotational_momentum_perstate = -hmm_1stlevel.cycle_metrics.rotational_momentum_perstate;
    
    % controls
    if whichstudy==1
        hmm_1stlevel.control.quintile.rotational_momentum = -hmm_1stlevel.control.quintile.rotational_momentum;
        hmm_1stlevel.control.ntile.rotational_momentum = -hmm_1stlevel.control.ntile.rotational_momentum;
        hmm_1stlevel.control.cardiac.rotational_momentum = -hmm_1stlevel.control.cardiac.rotational_momentum;
        
        for k=1:9
            hmm_1stlevel.control.timelagged.simulation{k}.cycle_metrics.rotational_momentum = -hmm_1stlevel.control.timelagged.simulation{k}.cycle_metrics.rotational_momentum;
            hmm_1stlevel.control.timelagged.simulation{k}.cycle_metrics.rotational_momentum_perstate = -hmm_1stlevel.control.timelagged.simulation{k}.cycle_metrics.rotational_momentum_perstate;
        end
        
    end
    
    % simulations
    for k=1:length(hmm_1stlevel.simulation)
        hmm_1stlevel.simulation{k}.cycle_metrics.rotational_momentum = -hmm_1stlevel.simulation{k}.cycle_metrics.rotational_momentum;
        hmm_1stlevel.simulation{k}.cycle_metrics.rotational_momentum_perstate = -hmm_1stlevel.simulation{k}.cycle_metrics.rotational_momentum_perstate;
    end
    hmm_1stlevel.simulation_average.cycle_metrics.rotational_momentum = -hmm_1stlevel.simulation_average.cycle_metrics.rotational_momentum;
    hmm_1stlevel.simulation_average.cycle_metrics.rotational_momentum_perstate = -hmm_1stlevel.simulation_average.cycle_metrics.rotational_momentum_perstate;
    
    
    % perm
    hmm_1stlevel.perm.rotational_momentum = -hmm_1stlevel.perm.rotational_momentum;
    
    
    %% redo
    cfg=[];
    cfg.method = 'montecarlo';
    cfg.statistic = 'depsamplesT';
    cfg.design = [ones(1,config.nSj), 2*ones(1,config.nSj); 1:config.nSj, 1:config.nSj];
    cfg.ivar = 1;
    cfg.uvar = 2;
    cfg.numrandomization = 100000;
    cfg.correcttail = 'prob';
    cfg.tail = 1;
    
    dat1=[];
    dat1.dimord = 'rpt_chan_time';
    dat1.label{1} = 'metric';
    
    measures = {'rotational_momentum', 'rotational_momentum_perstate'};
    for im = measures
        m = im{1};
        
        dat1.time=1:size(hmm_1stlevel.cycle_metrics.(m),2);
        dat1.trial = [];
        dat2=dat1;
        
        dat1.trial(:,1,:) = hmm_1stlevel.cycle_metrics.(m);
        dat2.trial(:,1,:) = hmm_1stlevel.simulation{1}.cycle_metrics.(m);
        
        hmm_1stlevel.metric_vs_sim.(m) = ft_timelockstatistics(cfg, dat1, dat2);
        save
        dat2.trial(:,1,:) = hmm_1stlevel.simulation_average.cycle_metrics.(m);
        hmm_1stlevel.metric_vs_sim_avg.(m) = ft_timelockstatistics(cfg, dat1, dat2);
    end
    save(config.metricfile, 'hmm_1stlevel', '-append')
    sprintf('Study %d Done', whichstudy)
    
end
