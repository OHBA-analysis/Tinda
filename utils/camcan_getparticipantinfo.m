function [subjinfo,subjid] = camcan_getparticipantinfo(config)

if config.nSj==600 % Matlab Fit
    participantinfo = tdfread(config.participantfile);

    basedir = '/ohba/pi/mwoolrich/datasets/CamCan_2021/';
    load(config.matfilelist);
    
    subjinfo = zeros(config.nSj,2);
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
    fname = fullfile(basedir, 'ParticipantCovariates/RTsimple_summary.xlsx');
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
    
    
    
elseif config.nSj==612 % Python fit
    q=tdfread('/ohba/pi/mwoolrich/mvanes/Projects/Tinda/Study4/ChetFit/camcan_all.csv');
    q=cellstr(q.x0x2CID0x2C0x22Sex_0x2810x3Dmale0x2C_20x3Dfemale0x290x220x2CAge);
    
    participantinfo = tdfread(fullfile('/ohba/pi/mwoolrich/datasets/CamCan_2021', 'ParticipantCovariates/participants.tsv'));
    
    cell_dat_split = cellfun(@(x)regexp(x,',','split'),q,'UniformOutput',0);
    
    for k=1:numel(cell_dat_split)
        subjinfo(k,:) = nan(1,5);
        a=cell_dat_split{k};
        subjinfo(k,1:3) = cellfun(@str2num, a([4,5,3]));
        subjid{k} = a{2};
        
        ix = find(~cellfun(@isempty,strfind(cellstr(participantinfo.participant_id), a{2})));        
        subjinfo(k,4) = participantinfo.tiv_cubicmm(ix);
    end
end