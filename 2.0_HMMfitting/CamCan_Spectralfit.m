%% and fit wavelet estimation:
load('/Volumes/CamsHD2/CamCan_2021/HMM/matfiles/filelist.mat')


%%
WTSpectfolder = '/Volumes/CamsHD2/CamCan_2021/HMM/WTspect/'
if ~isdir(WTSpectfolder)
    mkdir(WTSpectfolder)
end
hmmfile = ['/Volumes/CamsHD2/CamCan_2021/HMM/hmm_analysis1.mat'];
load(hmmfile,'hmm');

% set wavelet options:
options_wt = struct('Fs',250); % Sampling rate - for the 25subj it is 300
options_wt.fpass = [1 45];  % band of frequency you're interested in
options_wt.p = 0; %0.01; % interval of confidence  
options_wt.to_do = [1 0]; % turn off pdc
options_wt.order = 0;
options_wt.embeddedlags = -7:7;
d = length(options_wt.embeddedlags) - 1; 

for i=337:length(mat_files_orth)
    load(mat_files_orth{i});
    
    % get soft state timecourses:
    fprintf(['\nDecoding subj: ',int2str(i)]);
    Gamma = hmmdecode(X,T,hmm,0);
    save(mat_files_orth{i},'Gamma','-append');
    
    % and fit wavelet psd estimate:
    fprintf(['\n RUNNING WAVELET ON SUBJECT: ',int2str(i),'\n']);
    fitmt_subj = hmmspectrawt(X,T,Gamma,options_wt);
    fitmt_subj.state = rmfield(fitmt_subj.state,'ipsd');
    fitmt_subj.state = rmfield(fitmt_subj.state,'pcoh');
    fitmt_subj.state = rmfield(fitmt_subj.state,'phase');
    outputfilewt = [WTSpectfolder 'hmm_spectrawt_sj', int2str(i),'.mat'];
    save(outputfilewt,'fitmt_subj')
end
