function [info,labels] = MEGUK_getparticipantinfo()

behav_dir = '/Users/chiggins/data/Neuron2020/CanonicalRS/PsychoMet_and_Info/';
%load('/Users/chiggins/data/Neuron2020/CanonicalRS/PsychoMet_and_Info/PsychoMetsandID.txt')

M = dlmread([behav_dir,'PsychoMetsandID.txt'],' ',0,0);
fname_subs = [behav_dir,'RS_subs.mat'];
if ~isfile(fname_subs)
    % find the subjects that include eye_open:
    workingdir = '/Volumes/CamsHD2/NottinghamRS/raw_data/';
    datadir = workingdir;

    subjectdirs = dir([datadir '/3*']);
    subjectdirs = sort_nat({subjectdirs.name});

    for s=1:71
        dsfile = dir([datadir subjectdirs{s} '/*Eyes_Open*.ds']); 
        includerest(s) = ~isempty(dsfile);
    end
    validsubjs = find(includerest);
    for s=1:55
        goodsubs(s) = str2num(subjectdirs{validsubjs(s)});
    end

    save(fname_subs,'goodsubs');
else
    load(fname_subs);
end

% eliminate rows of M that have no RS:
info = [];
for i=1:length(goodsubs)
    r_ind = find(M(:,1) == goodsubs(i));
    if isempty(r_ind) || length(r_ind)>1
        error('Blank or duplicate entry')
    end
    info = [info;M(r_ind,:)];
end

labels = {'ID','IQ','Average Sternberg RT (across Loads)','Total Sternberg Accuracy',...
    'IPIP- Extraversion','IPIP- Agreeableness','IPIP- Conscientiousness','IPIP- Emotional Stability',...
    'IPIP- Intelligence','SPQ Score Full (based on P.Liddle Suggestion)','Age (years 2dec places, at time of MEG Scan)',...
    'Gender (1= Male 2=Female)','Open Rest Data Suitable for Group Analysis (1 = Yes, 0= No)',...
    'Closed Rest Data Suitable for Group Analysis (1 = Yes, 0= No)','SITI Data Suitable for Group Analysis (1 = Yes, 0= No)',...
    'LITI Data Suitable for Group Analysis (1 = Yes, 0= No)','Gamma Data Suitable for Group Analysis (1 = Yes, 0= No)',...
    'Sternberg Data Suitable for Group Analysis (1 = Yes, 0= No)','CEST Data Collected',...
    'CEST Data Suitable for Further Analysis','HADS: Depression','HADS: Anxiety','RRS: Rumination',...
    'SIPI: Daydreaming Full Score','Handedness (-12 = fully left, 0 = ambi, 12 = right)','fMRI Data Collected?'}

end