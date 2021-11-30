% exploratory poisson analysis of cycle timecourse

% model cycle timecourse as Poisson process:
clear all;
whichstudy = 4;
config = getStudyDetails(whichstudy);
load(config.metricfile)
color_scheme = set1_cols;
%% create timecourse of sequence emission per subject:
 
for subnum=1:config.nSj
    fprintf(['\nSubj',int2str(subnum)]);
    if whichstudy==4
        temp = hmm_2ndlevel.cyctimes_all{subnum}{1};
    else
        temp = hmm_2ndlevel.cyctimes_all{subnum};
    end
    
    cyctimes = cumsum(temp);
    timecourse = zeros(cyctimes(end),1);
    timecourse(cyctimes) = 1;
    cyctimes = cumsum(temp(randperm(length(temp))));
    timecourse_shuffle = zeros(cyctimes(end),1);
    timecourse_shuffle(cyctimes) = 1;
    for iW = 1:10*config.sample_rate
        count = conv(timecourse,ones(iW,1),'same');
        fanofactor(subnum,iW) = var(count) ./ mean(count);
        count = conv(timecourse_shuffle,ones(iW,1),'same');
        fanofactor_null(subnum,iW) = var(count) ./ mean(count);
    end
end
%%
figure();
shadedErrorBar((1:10*config.sample_rate)./config.sample_rate,mean(fanofactor),std(fanofactor)./sqrt(config.nSj),{'LineWidth',2});
plot4paper('Window size (seconds)','Fano Factor');
hold on;
plot((1:10*config.sample_rate)./config.sample_rate,mean(fanofactor_null),'k--');
print([config.figdir ,'supp/CycleRate_FanoFactor'],'-dpng');
%%
figure();
Sjs_to_plot = randi(config.nSj,1,12);
to_plot = hmm_2ndlevel.cyctimes_all(Sjs_to_plot);
for i=1:length(to_plot)
   to_plot{i} =  to_plot{i}{1}./config.sample_rate;
end
distributionPlot(to_plot,'showMM',2,'color',{color_scheme{1:12}});
ylim([0,10]);
print([config.figdir ,'supp/CycleRate_randomsubjects'],'-dpng');
%% save to file for later analysis:


% save metrics for later analysis:
if exist(config.metricfile)
    load(config.metricfile,'supp');
    supp.fanofactor = fanofactor;
    save(config.metricfile,'supp','-append')
else
    error('Cant find config file')
end
