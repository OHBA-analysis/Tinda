% this script is an alternative to the secondLevelHMMfit. Instead of
% running the 2nd level HMM on the vpath in order to get cycle statistics,
% this script replaces the vpath by the state's position on the unit
% circle, after which it runs spectral analysis on this. 

if ~exist('whichstudy','var')
    whichstudy = 1; % 1 denotes the hmm model run in Higgins2020_neuron
end
config = getStudyDetails(whichstudy);
% other preliminary setup for plotting etc:
color_scheme = set1_cols();

%% Set Poisson window length to average state lifetime:

temp = load(fullfile(config.hmmfolder,config.hmmfilename));

hmm = temp.hmm;
if ~isfield(hmm,'gamma') && whichstudy<4
    hmm.gamma = temp.Gamma;
    hmm.statepath = temp.vpath;
end
if whichstudy<3
    %load(config.prepdatafile,'hmmT','subj_inds');
    hmm = hmm_permutestates(hmm,temp.new_state_ordering);
    for i=1:config.nSj
        hmmT{i} = sum(hmm.subj_inds==i);
    end
elseif whichstudy==3 
    hmmT = temp.T_all;
    hmm.subj_inds = zeros(size(hmm.statepath));
    t_offset = 0;
    for i=1:length(hmmT)
        t_length = sum(hmmT{i}) - length(hmmT{i})*(length(hmm.train.embeddedlags)-1);
        hmm.subj_inds(t_offset + [1:t_length]) = i;
        t_offset = t_offset + t_length;
    end
    hmm.subj_inds = ceil(hmm.subj_inds/3); % account for multiple runs per subj
    hmmTold = reshape(hmmT,3,config.nSj);
    hmmTsubj = cell(config.nSj);
    for i=1:config.nSj
        hmmTsubj{i} = [hmmTold{1,i},hmmTold{2,i},hmmTold{3,i}]-(length(hmm.train.embeddedlags)-1);
    end
    hmmT = hmmTsubj;
    clear hmmTsubj hmmTold;
elseif whichstudy==4
    hmmT = temp.T_all;
    % correct for embedded lags:
    for i=1:length(hmmT)    
        hmmT{i} = hmmT{i} - (length(hmm.train.embeddedlags)-1);
    end
    load(config.matfilelist);
end
clear temp vpath;

if whichstudy<4
    Gamma = hmm.gamma;
end
K = hmm.K;
opts = [];
opts.K = 12;
opts.Fs = config.sample_rate;
FO = zeros(K,K,2,config.nSj);
for subnum=1:config.nSj
    if whichstudy<4
        vpath{subnum} = hmm.statepath(hmm.subj_inds==subnum);
    else
        temp = load(mat_files_orth{subnum},'vpath');
        vpath{subnum} = temp.vpath;
    end
    LT = getStateLifeTimes(vpath{subnum},hmmT{subnum},opts);
    LT_sj(subnum,:) = cellfun(@mean,LT);
end

Poisswindow = ceil(mean(LT_sj(:)));

%% 
W = Poisswindow; % window length to sum over
W = 125; %set arbitrary half second window
Noverlap = 1; % number of points to overlap over successive windows
T_poiss = [];
X_poiss = [];
Poiss_subj_inds = [];
t_last = 0;

% get the circle plot positions
optimalseqfile = [config.hmmfolder,'bestseq',int2str(whichstudy),'.mat'];
if ~isfile(optimalseqfile)
    bestsequencemetrics = optimiseSequentialPattern(FO);
    save(optimalseqfile,'bestsequencemetrics');
else
    load(optimalseqfile);
end
bestseq = bestsequencemetrics{2};
% save each subject's rotational strength by this metric:
disttoplot_manual = zeros(12,1);
for i3=1:12
    disttoplot_manual(bestseq(i3)) = exp(sqrt(-1)*i3/12*2*pi);
end

if Noverlap~=0
    overlapstring='_overlappingWindows';
end
for iSj=1:config.nSj
    %vpTemp = vpath(hmm.subj_inds==iSj,:);
    %sjind = 1+ (iSj-1)*3;
    %TSj{iSj} = cell2mat(T_all(sjind:(sjind+2)));
    %TSj{iSj} = TSj{iSj} - 14; % account for delay embedding
    vpTemp = vpath{iSj};%((t_last+1) : t_last+sum(TSj{iSj}));
    % replace the state count by its position on the unit circle (according
    % to the best sequence)
    vpcircle = vpTemp;
    for k=1:12
      vpcircle(vpcircle==k) = disttoplot_manual(k);
    end
    t_last = t_last + sum(hmmT{iSj});
    fprintf(['\nProcessing subject ',int2str(iSj)]);
    if whichstudy==4 %reset for each subject and save, rather than concatenate
        X_poiss = []; 
        T_poiss = [];
    end
    for i=1:length(hmmT{iSj})
        % Event count over window of length W:
        if Noverlap==0
            T_poiss = [T_poiss;floor(hmmT{iSj}(i)/W)];
            on_t = 1 + sum(hmmT{iSj}(1:(i-1)));
        else
            % set for overlapping windows
            T_poiss = [T_poiss;hmmT{iSj}(i)-floor(W/2)];
            on_t = sum(hmmT{iSj}(1:(i-1)));
        end
        
            
        if T_poiss(end)>1
            if Noverlap==0
                temp = reshape(vpTemp(on_t:(W*T_poiss(end)+on_t-1),:),[W,T_poiss(end)]);
            else
                temp = zeros(W,T_poiss(end));
                Xcircle = temp;
                for i2=1:(T_poiss(end)-floor(W/2))
                    temp(:,i2) = vpTemp(on_t + [i2:(i2+W-1)]);
                    Xcircle(:,i2) = vpcircle(on_t + [i2:(i2+W-1)]);
                end
            end
            temp2 = zeros(T_poiss(end),hmm.K);
            for k=1:hmm.K
                temp2(:,k)=sum(temp==k,1);
            end
            Xcircle = angle(sum(Xcircle,1));
            X_poiss = [X_poiss;temp2];
            Poiss_subj_inds = [Poiss_subj_inds;repmat(iSj,T_poiss(end),1)];
        else
            T_poiss(end)=[];
        end
    end
    T_poiss(T_poiss==0)=[];
    if whichstudy==4
        if ~isfolder([config.hmmfolder,'Poissdata_',int2str(W),overlapstring])
            mkdir([config.hmmfolder,'Poissdata_',int2str(W),overlapstring]);
        end
        mat_files_poiss{iSj} = [config.hmmfolder,'Poissdata_',int2str(W),overlapstring,'/PoissDataSj',int2str(iSj)];
        X = X_poiss;
        T = T_poiss;
        save(mat_files_poiss{iSj},'X','T');
    end
end


if whichstudy==4
    save([config.hmmfolder,'Poissdata_',int2str(W),overlapstring,'/filelist.mat'],'mat_files_poiss');
end
% % optionally for large datasets save to individual files:
% n_runs = 10;
% if whichstudy==4
%     % need to save to individual matfiles for memory constraints
%     mkdir([config.hmmfolder,'Poissdata_',int2str(W)]);
%     t0=1;
%     for subnum=1:config.nSj
%         fprintf(['\nSaving Sj ',in2str(subnum)]);
%         X = X_poiss(Poiss_subj_inds==subnum,:);
%         if subnum<config.nSj
%             T = T_poiss(t0:(find(cumsum(T_poiss)==find(Poiss_subj_inds==subnum+1,1)-1)));
%             t0 = t0+length(T);
%         else
%             T = T_poiss(t0:end);
%         end
%         mat_files_poiss{subnum} = [config.hmmfolder,'Poissdata_',int2str(W),'/PoissDataSj',int2str(subnum)];
%         save(mat_files_poiss{subnum},'X','T');
%         T_cells{subnum} = T;
%     end
%     n_runs=1;
%     clear X temp temp2 vpath vpTemp X_poiss %T_poiss % save memory
% end

%%
n_runs = 1;
for i_run = 1:n_runs
options = [];
options.K = 3; % Note this didn't work with 4 states; one is knocked out and results in wierd behaviour
options.distribution = 'poisson';
options.Pstructure = eye(options.K) + diag(ones(1,options.K-1),1);
options.Pstructure(options.K,1) = 1;
options.initrep = 4; % this the number of parallel cores

% NOTE: initialisation is shite - winds up with states that are extremely
% rarely on. just do random init here:

%options.Gamma = rand(T_poiss(iSj),options.K);
%options.DirichletDiag = sum(T_poiss);
Sjs_init = randperm(config.nSj,30);
if whichstudy<4
    X_init = []; T_init = [];
    for iSj = Sjs_init
        X_init = [X_init;X_poiss((1+sum(T_poiss(1:(iSj-1)))):sum(T_poiss(1:iSj)),:)];
        T_init = [T_init;T_poiss(iSj)];
    end
    [hmminit,gammainit] = hmmmar(X_init,T_init,options);
%     fe_lowest = Inf;
%     for i_init=1:10
%         options.decodeGamma = 0;
%         [hmminit_temp,~,~,~,~,~,fe] = hmmmar(X_poiss,T_poiss,options);
%         if i_init==1 | fe(end)<fe_lowest
%             hmminit = hmminit_temp;
%             fe_lowest = fe(end);
%         end
%     end
else
    
    X_poiss = []; T_poiss = [];
    for iSj = Sjs_init
        load(mat_files_poiss{iSj});
        X_poiss = [X_poiss;X];
        T_poiss = [T_poiss;T];
    end
    fe_lowest = Inf;
    for i_init=1:10
        options.decodeGamma = 0;
        [hmminit_temp,~,~,~,~,~,fe] = hmmmar(X_poiss,T_poiss,options);
        if i_init==1 | fe(end)<fe_lowest
            hmminit = hmminit_temp;
            fe_lowest = fe(end);
        end
    end
end


options = rmfield(options,'initrep')
options.hmm = hmminit;

% % or set to figure 8 pattern:
% options.Pstructure(1,4)=1;
% options.Pstructure(3,1)=1;
% options.Pstructure(3,4)=0;
options.decodeGamma = false;
options.standardise = false;
if whichstudy~=4
    [hmmtemp,Gammatemp,~,~,~,~,fehist] = hmmmar(X_poiss,T_poiss,options);

    if i_run==1 || fehist(end)<lowestfe
        hmmPoiss = hmmtemp;
        GammaPoiss = Gammatemp;
        lowestfe = fehist(end);
    end
    % record a few stats to get a sense of how deviant these really are:
    feall(i_run,1) = fehist(1);
    feall(i_run,2) = fehist(end);
    gamsum(i_run,:) = mean(Gammatemp);
else
    % set stochastic inference options:
%     options.BIGNinitbatch = 5;
%     options.BIGNbatch = 5;
%     options.BIGtol = 1e-7;
%     options.BIGcyc = 500;
%     options.BIGundertol_tostop = 5;
%     options.BIGdelay = 5;
%     options.BIGforgetrate = 0.7;
%     options.BIGbase_weights = 0.9;
    
    
    options.cyc=10; % limit the number of inference cycles to run on each noisy update
    %[hmmtemp] = hmmmar(mat_files_poiss,T_cells,options);
    % Hack to implemnet stochastic inference:
    n_batch = 30;
    for i=1:100
        thisbatch = sort(randperm(600,n_batch));
        X_poiss = []; T_poiss = [];
        % loading batch files:
        fprintf(['\nLoading batch files, batch',int2str(i)]);
        for j=1:n_batch
            temp = load(mat_files_poiss{thisbatch(j)});
            X_poiss = [X_poiss;temp.X];
            T_poiss = [T_poiss;temp.T];
        end
        fprintf(['\nFitting hmm to new batch:']);
        [hmmtemp,Gammatemp,~,~,~,~,fehist] = hmmmar(X_poiss,T_poiss,options);
        
        % update key params:
        alpha = 5;beta = 0.7;
        LR = 1-(i+alpha).^-beta;
        for k=1:3
            options.hmm.state(k).W.W_shape = LR*options.hmm.state(k).W.W_shape + (1-LR)*hmmtemp.state(k).W.W_shape;
            options.hmm.state(k).W.W_mean = LR*options.hmm.state(k).W.W_mean + (1-LR)*hmmtemp.state(k).W.W_mean;
            options.hmm.state(k).W.W_rate = LR*options.hmm.state(k).W.W_rate + (1-LR)*hmmtemp.state(k).W.W_rate;
        end
        options.hmm.Pi = LR*options.hmm.Pi + (1-LR)*hmmtemp.Pi;
        options.hmm.Pe = LR*options.hmm.Pe + (1-LR)*hmmtemp.Pe;
        options.hmm.P = LR*options.hmm.P + (1-LR)*hmmtemp.P;
        options.hmm.Dir_alpha = LR*options.hmm.Dir_alpha + (1-LR)*hmmtemp.Dir_alpha;
        options.hmm.Dir2d_alpha = LR*options.hmm.Dir2d_alpha + (1-LR)*hmmtemp.Dir2d_alpha;
        
        % monitor some params to see how they evolve:
        Wrec(i,:) = [options.hmm.state(k).W.W_mean,options.hmm.state(k).W.W_mean,options.hmm.state(k).W.W_mean];
    end
end
end
%%
if whichstudy~=4
    hmmPoiss.gamma = GammaPoiss;
    disttoplot = plotMDS_states(hmmPoiss);

    save([config.hmmfolder,'secondLevelHMM_Poiss_window',num2str(W),'_K',int2str(options.K),overlapstring,'.mat'],'hmmPoiss','feall','GammaPoiss','T_poiss','Poiss_subj_inds');
    % [~,new_state_orderingPoiss] = sort(disttoplot(:,1));
    % hmmPoiss = hmm_permutestates(hmmPoiss,new_state_orderingPoiss);
    % GammaPoiss = hmmPoiss.gamma;

    % plot convergence details;
    figure();subplot(3,1,1);
    plot(feall,'LineWidth',2);
    plot4paper('Init run','Free energy');
    legend('On init','On convergence')
    subplot(3,1,2);
    plot(std(gamsum,[],2),'*');
    plot4paper('Init run','Gamma std');
    subplot(3,1,3);
    plot(min(gamsum,[],2),'*');
    plot4paper('Init run','Min Gamma');
    %print([figdir,'Fig0_ConvergenceRecord_window',int2str(iW)],'-dpng');
else
    hmmPoiss = rmfield(options.hmm,'Gamma');
    save([config.hmmfolder,'secondLevelHMM_stoch_Poiss_window',num2str(W),'_K',int2str(options.K),overlapstring,'.mat'],'hmmPoiss');
    
    % and infer each subject's associated state timecourse:
    options.updateObs = 0;
    options.decodeGamma = 1;
    options.cyc=1;
    options.hmm = hmmPoiss;
        
    for i=1:600
        fprintf(['\n inferring STC for sj ',int2str(i)]);
        temp = load(mat_files_poiss{i});
        [~,Gamma] = hmmmar(temp.X,temp.T,options);
        save(mat_files_poiss{i},'Gamma','-append');
    end
end

%% finally, save summary stats to file for later analysis:
W=125;K=3;
overlapstring='_overlappingWindows';

% infer cycle times:
if whichstudy<4
    load([config.hmmfolder,'secondLevelHMM_Poiss_window',num2str(W),'_K',int2str(K),overlapstring,'.mat'],'hmmPoiss','feall','GammaPoiss','T_poiss','Poiss_subj_inds');

    samp_minute = (config.sample_rate*60); % split into minute by minute chunks
    for subnum = 1:config.nSj

        Gamtemp = hmmPoiss.gamma(Poiss_subj_inds==subnum,:);
        if whichstudy==3
            cycletimes = getStateCycleTimes(Gamtemp,T_poiss([1:3]+(subnum-1)*3));
            lifetimes_meta{subnum} = getStateLifeTimes(Gamtemp,T_poiss([1:3]+(subnum-1)*3));
        else
            cycletimes = getStateCycleTimes(Gamtemp,T_poiss(subnum));
            lifetimes_meta{subnum} = getStateLifeTimes(Gamtemp,T_poiss(subnum));
        end
        cycletime_mu(subnum,:) = mean(cycletimes);
        cycletime_std(subnum,:) = std(cycletimes);
        cycletime_med(subnum,:) = median(cycletimes);
        FO_meta(subnum,:) = getFractionalOccupancy(Gamtemp,length(Gamtemp));
        cyctimes{subnum} = cycletimes;
        
        if whichstudy==3 % get session by session detail:
            for subnum = 1:config.nSj
                for isession=1:3
                    t0 = 1+sum(T_poiss(1:((subnum-1)*3)+(isession-1)));tend = sum(T_poiss(1:((subnum-1)*3)+(isession)));
                    Gamtemp = hmmPoiss.gamma(t0:tend,:);
                    cycletimes = getStateCycleTimes(Gamtemp,tend-t0+1);
                    cycletime_mu_sess(subnum,isession) = mean(cycletimes);
                    cycletime_std_sess(subnum,isession) = std(cycletimes);
                    cycletime_med_sess(subnum,isession) = median(cycletimes);
                    FO_meta_sess(subnum,isession,:) = getFractionalOccupancy(Gamtemp,length(Gamtemp));
                end
            end
        end
%         % also get minute by minute detail:
%         for imin=1:5
%             startseg = find(cumsum(T)>(imin-1)*samp_minute,1)-1;
%             endseg = find(cumsum(T)>(imin)*samp_minute,1)-1;
%             if ~isempty(startseg) && ~isempty(endseg)
%                 if startseg==endseg
%                     T_sub = samp_minute;
%                 else
%                     T_sub = sum(T(1:startseg+1)) - (imin-1)*samp_minute;
%                     T_sub = [T_sub;T(startseg+2:endseg)];
%                     T_sub = [T_sub;samp_minute - sum(T_sub)];
%                 end
%                 try
%                     Gamtemp = Gamma((imin-1)*samp_minute + [1:samp_minute],:);
%                 catch
%                     Gamtemp = Gamma((imin-1)*samp_minute :end,:);
%                     T_sub = T_sub(cumsum(T_sub)<length(Gamtemp));
%                     T_sub = [T_sub;length(Gamtemp)-sum(T_sub)];
%                 end
%                 temp = getStateCycleTimes(Gamtemp,T_sub);
%                 cyctime{subnum,imin} = temp;
%                 cycletime_mu_min(subnum,imin) = mean(temp);
%                 cycletime_med_min(subnum,imin) = median(temp);
%                 cycletime_std_min(subnum,imin) = std(temp);
%             end
%         end 
    end
else
    load([config.hmmfolder,'secondLevelHMM_stoch_Poiss_window',num2str(W),'_K',int2str(K),overlapstring,'.mat'],'hmmPoiss','feall','GammaPoiss','T_poiss','Poiss_subj_inds');

    if ~contains(config.Poiss_dir,'overlappingWindows')
        figdir = [config.figdir,'4_covariates_W',int2str(W),'/'];
        mkdir(figdir);
        clear FO 
        for subnum=1:600
            Gamtemp = hmmPoiss.gamma(Poiss_subj_inds==subnum,:);
            cycletimes = getStateIntervalTimes(Gamtemp,length(Gamtemp));
            cycletime_mu(subnum,:) = cellfun(@mean,cycletimes);
            cycletime_std(subnum,:) = cellfun(@std,cycletimes);
            cycletime_med(subnum,:) = cellfun(@median,cycletimes);
            FO(subnum,:) = getFractionalOccupancy(Gamtemp,length(Gamtemp));
            cyctimes{subnum} = cycletimes;
            lifetimes_meta{subnum} = getStateLifeTimes(Gamtemp,length(Gamtemp));
        end

    else
        figdir = [config.figdir,'4_covariates_W',int2str(W),'_overlappingWindows/'];
        mkdir(figdir);
        clear cycletimes cycletime_mu cycletime_std cycletime_med FO cyctimes lifetimes cycletime_mu_min cycletime_med_min cycletime_std_min
        load([config.Poiss_dir,'filelist.mat'])
        samp_2minute = config.sample_rate*2*60;
        for subnum=1:600
            fprintf(['\nSubj: ',int2str(subnum)]);
            load(mat_files_poiss{subnum},'Gamma','T');
            %cycletimes = getStateIntervalTimes(Gamma,T);
            cycletimes{1} = getStateCycleTimes(Gamma,T);
            cycletime_mu(subnum,:) = cellfun(@mean,cycletimes);
            cycletime_std(subnum,:) = cellfun(@std,cycletimes);
            cycletime_med(subnum,:) = cellfun(@median,cycletimes);
            FO_meta(subnum,:) = getFractionalOccupancy(Gamma,length(Gamma));
            cyctimes{subnum} = cycletimes;
            lifetimes_meta{subnum} = getStateLifeTimes(Gamma,T);
            % also get minute by minute detail:
            for imin=1:5
                startseg = find(cumsum(T)>(imin-1)*samp_2minute,1)-1;
                endseg = find(cumsum(T)>(imin)*samp_2minute,1)-1;
                if ~isempty(startseg) && ~isempty(endseg)
                    if startseg==endseg
                        T_sub = samp_2minute;
                    else
                        T_sub = sum(T(1:startseg+1)) - (imin-1)*samp_2minute;
                        T_sub = [T_sub;T(startseg+2:endseg)];
                        T_sub = [T_sub;samp_2minute - sum(T_sub)];
                    end
                    try
                        Gamtemp = Gamma((imin-1)*samp_2minute + [1:samp_2minute],:);
                    catch
                        Gamtemp = Gamma((imin-1)*samp_2minute :end,:);
                        T_sub = T_sub(cumsum(T_sub)<length(Gamtemp));
                        T_sub = [T_sub;length(Gamtemp)-sum(T_sub)];
                    end
                    temp = getStateCycleTimes(Gamtemp,T_sub);
                    cyctime{subnum,imin} = temp;
                    cycletime_mu_min(subnum,imin) = mean(temp);
                    cycletime_med_min(subnum,imin) = median(temp);
                    cycletime_std_min(subnum,imin) = std(temp);
                end
            end 
        end
    end
end

hmm_2ndlevel = [];
hmm_2ndlevel.cyctime_mu = cycletime_mu;
hmm_2ndlevel.cyctime_med = cycletime_med;
hmm_2ndlevel.cyctime_std = cycletime_std;
hmm_2ndlevel.FO = FO_meta;
hmm_2ndlevel.cyctimes_all = cyctimes;
hmm_2ndlevel.LT_all = lifetimes_meta;

if whichstudy==3
    hmm_2ndlevel.cycletime_mu_sess = cycletime_mu_sess;
    hmm_2ndlevel.cycletime_std_sess = cycletime_std_sess;
    hmm_2ndlevel.cycletime_med_sess = cycletime_med_sess;
    hmm_2ndlevel.FO_meta_sess = FO_meta_sess;
elseif whichstudy==4
    hmm_2ndlevel.cyctime_min = cyctime;
    hmm_2ndlevel.cycletime_mu_min=cycletime_mu_min;
    hmm_2ndlevel.cycletime_med_min=cycletime_med_min;
    hmm_2ndlevel.cycletime_std_min=cycletime_std_min;
end
% save metrics for later analysis:
if exist(config.metricfile)
    save(config.metricfile,'hmm_2ndlevel','-append');
else
    save(config.metricfile,'hmm_2ndlevel')
end