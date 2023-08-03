%% CamCan sanity checks:
if isfolder('/Volumes/T5_OHBA/')
  basedir = '/Volumes/T5_OHBA/Projects/Tinda';
else
  basedir = '/ohba/pi/mwoolrich/datasets/CamCan2021/';
end
participantinfo = tdfread(fullfile(basedir, 'ParticipantCovariates/participants.tsv'));
camcandir_matfiles = fullfile(basedir, 'HMM/matfiles/');
load([camcandir_matfiles,'filelist'],'mat_files_orth','T_all');
% MvE addition:
for k=1:length(mat_files_orth)
  mat_files_orth{k} = replace(mat_files_orth{k}, [basedir, '/HMM/matfiles/'], camcandir_matfiles);
end
for i=1:length(mat_files_orth)
%     mat_files_orth{i}(1)='E';
    load(mat_files_orth{i},'vpath')
    for k=1:12
        FO(i,k) = mean(vpath==k);
        stONtimes = find(diff([0;vpath==k;0])==1);
        stOFFtimes = find(diff([0;vpath==k;0])==-1);
        LT(i,k) = mean(stOFFtimes-stONtimes);
        IT(i,k) = mean(stONtimes(2:end)-stOFFtimes(1:end-1));
    end
end

%% check for structure with age:
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
    end
end

%% see which states correlate with age:

for k=1:12
    [A,B] = corrcoef(FO(:,k),subjinfo(:,1));
    rho(k) = A(1,2);
    pval(k) = B(1,2);
end
pval

%% see whether mean LTs or ITs correlate with age:

figure();
subplot(1,2,1);
scatter(subjinfo(:,1),median(LT,2)/250);
[A,B] = corrcoef(mean(LT,2),subjinfo(:,1));
plot4paper('Age','State Lifetimes (sec)');
title(['rho=', num2str(A(2,1)),', p=',num2str(B(2,1))]);

subplot(1,2,2);
scatter(subjinfo(:,1),median(IT,2)/250);
[A,B] = corrcoef(mean(IT,2),subjinfo(:,1));
plot4paper('Age','State interval times (sec)');
title(['rho=', num2str(A(2,1)),', p=',num2str(B(2,1))]);

%% independent of HMM: fit pwelch to data and plot as function fo age:
clear psdrec
for i=1:600
    fprintf(['\nRunning on subj:',int2str(i)])
    load(mat_files_orth{i},'X');
    [temp,F] = pwelch(X,250,0,250,250);
    psdrec(i,:) = mean(temp,2);
end

%%

B = pinv([ones(600,1),normalise(subjinfo(:,1))])*psdrec;
figure();
subplot(1,2,1);
plot(F,B(1,:),'LineWidth',2);
grid on;xlim([0,40]);
plot4paper('Frequency (Hz)','Mean psd over subjects');
subplot(1,2,2);
plot(B(2,:),'LineWidth',2);
grid on;
xlim([0,40]);
plot4paper('Frequency (Hz)','Modulatory effect of age');
% compute significance:
for i=1:length(F)
    [~,~,~,~,stats] = regress(demean(psdrec(:,i)),[ones(600,1),normalise(subjinfo(:,1))]);
    pvals_psd(i) = stats(3);
end
Bsig = B(2,:);
Bsig(pvals_psd>0.05/40) = NaN;
hold on;
plot(Bsig,'LineWidth',4,'Color','Black');