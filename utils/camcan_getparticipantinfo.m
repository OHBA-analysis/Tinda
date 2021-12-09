function [subjinfo,subjid] = camcan_getparticipantinfo(config)
load(config.matfilelist);
participantinfo = tdfread(config.participantfile);
subjinfo = zeros(600,2);
for i=1:653
    temp = strfind(mat_files_orth,participantinfo.participant_id(i,:));
    %temp = cell2mat(temp);
    rowind = find(~cellfun(@isempty,temp));
    
    if ~isempty(rowind)
        subjinfo(rowind,1) = participantinfo.age(i);
        subjinfo(rowind,2) = participantinfo.hand(i);
        subjinfo(rowind,3) = participantinfo.gender_code(i);
        subjinfo(rowind,4) = participantinfo.tiv_cubicmm(i);
        subjid{rowind} = participantinfo.participant_id(i,:);
        
    end
    
end
% get RTs:
fname = '/ohba/pi/mwoolrich/datasets/CamCan_2021/ParticipantCovariates/RTsimple_summary.xlsx'
headerlines = 0;
DATA = importdata(fname,' ',headerlines);
for i=1:600
    thissub = find(contains(DATA.textdata(:,1),subjid{i}(7:end)));
    if isempty(thissub)
        RTmean(i,1) = NaN;
    else
        RTmean(i,1) = DATA.data(thissub-1,8);
    end
end
subjinfo(:,5) =  RTmean;
end
