whichstudy=1;
config=getStudyDetails(whichstudy)
load(fullfile(config.hmmfolder,config.hmmfilename));

q=load('/ohba/pi/mwoolrich/mvanes/Projects/Tinda/Study1/secondLevelHMM_Poiss_window125_K3_overlappingWindows.mat');

[~, statepath] = max(q.hmmPoiss.gamma,[],2);

% add in the window length so that there is correspondence with the first
% level
tmp = zeros(size(hmm.subj_inds));
cnt = 1;
for iSj=1:config.nSj
    ix = find(hmm.subj_inds==iSj);
    ix = ix(1:end-124);
    tmp(1,ix) = statepath(cnt:cnt+length(ix)-1);
    cnt = cnt+length(ix);
end
statepath = tmp;

hmm_copy = hmm;
hmm_copy.statepath = statepath;
