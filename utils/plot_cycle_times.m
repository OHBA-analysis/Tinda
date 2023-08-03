% script called to load viterbi paths inferred and hmm objects and run
% post-hoc sequence analysis:
if ~exist('whichstudy','var')
    whichstudy = 1; % 1 denotes the hmm model run in Higgins2020_neuron
end
config = getStudyDetails(whichstudy);
% other preliminary setup for plotting etc:
color_scheme = colorscheme(whichstudy);
hotcold = cmap_hotcold;
clr = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330]};

% include option for simulations here - results to be saved in adjacent
% folder:
simtests = false;
if simtests
    config.figdir = strrep(config.figdir,['Study',int2str(whichstudy)],['Study',int2str(whichstudy),'_simtest']);
    mkdir(config.figdir);
end
use_WB_nnmf=true; % whether or not to use the wide band NNMF Diego's Nature Comms to select power and coherence (alternative is selecting 1-30 Hz)
useMT = true; % use the MT results instead of Cam's wavelet approach



if strcmp(config.reordering_states, 'coherence')
    optimalseqfile = [config.resultsdir,'bestseq',int2str(whichstudy),'_coherence' ,'.mat'];
elseif strcmp(config.reordering_states, 'study1matched')
    optimalseqfile = [config.resultsdir,'bestseq',int2str(whichstudy),'_study1matched' ,'.mat'];
else
    optimalseqfile = [config.hmmfolder,'bestseq',int2str(whichstudy),'.mat'];
end
load(optimalseqfile);
bestseq = bestsequencemetrics{1};

load(config.metricfile)

color_scheme2 = {[179 225 153]/255, [140 132 126]/255, [236 123 119]/255, [245 191 134]/255};
color_scheme3 = {[0 0 0]};
%% get a figure with all timings of 1st and 2nd level
%% Get figure showing number of visits per metastate
% load 2ndlevel vpath
W=17;
K2=4;
overlapstring='_overlappingWindows';
load([config.resultsdir,'secondLevelHMM_Poiss_window',num2str(W),'_K',int2str(K2),overlapstring,'.mat'])

for subnum=1:config.nSj
  Gamtemp = hmmPoiss.gamma(Poiss_subj_inds==subnum,:);
        vpTemp = vpath{subnum};
        if length(T_poiss)~=config.nSj
            error('Need to realign state timecourses for sub segments')
        else
            vpTemp = vpTemp(ceil(W/2):end-floor(W/2));
            if length(vpTemp)~=length(Gamtemp)
                error('Sizes misaligned')
            end
        end

        if whichstudy==3
            num_transitions = getMicroStateCycleVisits(Gamtemp,T_poiss([1:3]+(subnum-1)*3),vpTemp);
        else
            [num_transitions,microsequences{subnum}] = getMicroStateCycleVisits(Gamtemp,T_poiss(subnum),vpTemp);
        end
        microtransitions_all{subnum} = num_transitions;
        microtransitions_mu(subnum,:) = mean(num_transitions);
        microtransitions_median(subnum,:) = median(num_transitions);
end
%%
disttoplot_manual = zeros(12,2);
for i=1:12
    temp = exp(sqrt(-1)*(i+2)/12*2*pi);
    disttoplot_manual(bestseq(i),:) = [real(temp),imag(temp)];
end
    %circleposition = exp(j*(pi/2-[0:11]*2*pi/12));
    %circleposition = circleposition(bestseq);
    circleposition = disttoplot_manual(:,1) + sqrt(-1)*disttoplot_manual(:,2);
phaseshiftcounter = zeros(config.nSj,12);
phases = linspace(-pi,5*pi/6,12);
jointcounter = zeros(12,12,config.nSj);
for subnum=1:config.nSj

    if whichstudy==4
        clear vpath;
        vpTemp = load(mat_files_orth{subnum},'vpath','T');
        vpath{subnum} = vpTemp.vpath;
    end
    transitions = find(diff(vpath{subnum})~=0);
    phaseshift = zeros(length(transitions),1);
    for i=1:length(transitions)
        phase_t = angle(circleposition(vpath{subnum}(transitions(i))));
        phase_tplus1 = angle(circleposition(vpath{subnum}(transitions(i)+1)));
        phaseshift(i) = phase_tplus1 - phase_t;
    end
    % unwrap phases onto same scale:
    phaseshift(phaseshift > pi-0.01) = phaseshift(phaseshift > pi-0.01)-2*pi;
    phaseshift(phaseshift < -pi-0.01) = phaseshift(phaseshift < -pi-0.01)+2*pi;
    phaseshiftcounter(subnum,:) = hist(phaseshift,phases);

    % track phase shifts for consecutive transitions:
    for i=1:12
        pointsin = phaseshift(1:end-1) >= phases(i) - 0.1 & phaseshift(1:end-1) <= phases(i) + 0.1;
        pointsin = find(pointsin);
        for i2 = 1:12
            jointcounter(i,i2,subnum) = sum(phaseshift(pointsin + 1) >= phases(i2) - 0.1 & phaseshift(pointsin + 1) <= phases(i2) + 0.1);
        end
    end
end
phaseshiftcounter = phaseshiftcounter ./ repmat(sum(phaseshiftcounter,2),1,12);
%%
fig =setup_figure([],2,.5);
ax(1) = axes('Position',[0.125 0.175 0.25 0.75])
distributionPlot([hmm_1stlevel.LT_med./250, hmm_2ndlevel.LT_mu./250 hmm_2ndlevel.cyctime_mu], 'showMM',2,'color', [color_scheme, color_scheme2, color_scheme3])
view([90 90])
for k=1:17
  if k<=12
    lb{k} = sprintf('Microstate %d', k);
  elseif k<17
    lb{k} = sprintf('Meta state %d', k-12);
  else
     lb{k} = 'Cycle time';
  end
end
xticks(1:17)
xticklabels(lb)
% set(gca,'XTickLabelRotation',45)
vline(12.5, '--k')
vline(16.5, '--k')


ylabel('Mean Life/Cycle times (s)')
% subplot(1,2,2)
ax(2) = axes('Position',[0.45 0.175 0.2 0.75])

distributionPlot(microtransitions_mu, 'showMM',2, 'color', clr(1:2))
ylabel({'Mean number of state visits per cycle'})
% set(gca,'XTickLabel',{'Total state \\newline\\ visits','Unique state visits'})
set(gca,'XTickLabel',{sprintf('%s\\newline%s\\newline%s','Total', 'state', 'visits'),sprintf('%s\\newline%s\\newline%s','Unique', 'state', 'visits')})
set(gca,'XTickLabelRotation',0)


  
% figure();
ax(3) = axes('Position',[0.725 0.175 0.25 0.75])

bar([phases pi],fliplr([(mean(phaseshiftcounter)) 0]))
box off
xlabel('Phase shift')
ylabel({'Proportion of microstate transitions'})
% plot4paper('Phase Shift',)
labs = {'-\pi','-5\pi/6','-2\pi/3','-\pi/2','-\pi/3','-\pi/6','0','\pi/6','\pi/3','\pi/2','2\pi/3','5\pi/6'}
set(gca,'XTick',phases,'XTickLabel',labs)
% title('Histogram of Microstate Transition Phase Shifts')

save_figure([config.figdir, 'figure_supp_hmm_stats/fig_life_cycle_times'])




