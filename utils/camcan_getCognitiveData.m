function [cogdata,labels] = camcan_getCognitiveData(config)
if config.nSj==600 % Matlab Fit
    
    load('/ohba/pi/mwoolrich/datasets/CamCan_2021/ParticipantCovariates/CogDatAll.mat');
    load(config.matfilelist);
    cogdata = zeros(600,16);
    for i=1:size(CogDatAll,1)
        temp = strfind(mat_files_orth,num2str(CogDatAll(i,16)));
        %temp = cell2mat(temp);
        rowind = find(~cellfun(@isempty,temp));
        if ~isempty(rowind)
            cogdata(rowind,:) = CogDatAll(i,:);
        end
    end
elseif config.nSj==612 % Python Fit
    q=tdfread([config.resultsdir, 'camcan_all.csv']);
    q=cellstr(q.x0x2CID0x2C0x22Sex_0x2810x3Dmale0x2C_20x3Dfemale0x290x220x2CAge);
    
    cell_dat_split = cellfun(@(x)regexp(x,',','split'),q,'UniformOutput',0);
    cogdata = zeros(config.nSj, 16);
    for k=1:numel(cell_dat_split)
        a=cell_dat_split{k};
        subjid = a{2};
        cogdata(k,:) = [cellfun(@str2num, a([3,4,6:end])), str2num(subjid(7:end))];
    end
end
% now labels; note the following:
% 1) FldIn   : fluid intelligence           f   ex
% 2) FacReg  : face recognition             f   emo
% 3) EmoRec  : emotional recognition        f   emo
% 4) MltTs   : multitask                    f   ex      (negatively correlated with others - indicates reaction times)
% 5) PicName : picture naming               f   lang
% 6) ProV    : proverb comprehension        c   lang
% 7) MRSp    : motor speed                  f   mot     (negatively correlated with others - indicates reaction times)
% 8) MRCv    : motor speed varianve         f   mot
% 9) SntRec  : sentence comprehension       c/f lang
% 10) VSTM   : visual short-term memory     f   mem
% 11) StrRec : story recall                 f   mem
% 12) StW    : spot the word                c   lang
% 13) VrbFl  : verbal fluency               f   lang
% col1: M/F; col2: age; cols3-15: cognitive measures; col16: CCID
labels = {'Sex','Age','FldIn','FacReg','EmoRec','MltTs','PicName','ProV','MRSp','MRCv','SntRec',...
    'VSTM','StrRec','StW','VrbFl','CCID'};