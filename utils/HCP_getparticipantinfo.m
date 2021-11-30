function [info,headers] = HCP_getparticipantinfo(config)
    subjdata = readtable('/Users/chiggins/data/HCPAnalysis/behav/unrestricted_aquinn501_4_7_2017_9_4_13.csv');
    temp = readtable('/Users/chiggins/data/HCPAnalysis/behav/MEGfnames.csv');
    subj_ids = [];
    for i=1:size(temp,1)
        subj_id = str2num(temp{i,1}{1}(7:12));
        if ~isempty(subj_id) && length(intersect(subj_id,subj_ids))==0
            subj_ids = [subj_ids;subj_id];
        end
    end

    inds = [];
    for i=1:length(subj_ids)
        inds(i) = find(subjdata.Subject == subj_ids(i));
    end
    subjdata = subjdata(inds,:);

    % also load more detailed data and align:
    subjdata_detailed = readtable('/Users/chiggins/data/HCPAnalysis/behav/vars.txt');
    %headers = readtable('/Users/chiggins/data/HCPAnalysis/behav/column_headers.txt');
    clear headers;
    fid = fopen('/Users/chiggins/data/HCPAnalysis/behav/column_headers.txt');
    tline = fgetl(fid);
    i=1;
    while ischar(tline)
        temp = strrep(tline,' ','');
        temp = strrep(temp,'-','');
        headers{i,1} = temp;
        tline = fgetl(fid);i=i+1;
    end
    fclose(fid);

    inds = [];
    for i=1:length(subj_ids)
        inds(i) = find(subjdata_detailed.Var1 == subj_ids(i));
    end
    subjdata_detailed = subjdata_detailed(inds,:);
    subjdata_detailed.Properties.VariableNames = headers;
    
    info = table2array(subjdata_detailed);
    
    reactiontimevars = contains(headers,'RT');
    RTs = info(:,reactiontimevars);
    RTs = RTs(:,~any(isnan(RTs)));
    [~, RT_pc, eigspect]= pca(normalise(RTs));
    info = [info,RT_pc];
end