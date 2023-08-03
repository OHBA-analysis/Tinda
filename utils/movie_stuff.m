% misc script: take viterbi path, and compute fractional occupancy during
% periods between state 1 and 12 (directional)

basedir = '/ohba/pi/mwoolrich/datasets/ReplayData/';
wd = basedir;%[basedir,'data/Neuron2020/'];
% Add netlab and fmt
addpath( fullfile(osldir,'ohba-external','netlab3.3','netlab') );
addpath( fullfile(osldir,'ohba-external','fmt') );

download_path = [basedir,'Documents/MATLAB/Neuron2020/'];
addpath(genpath(download_path));

TESTRUN = false; % when this is on, just run full script for one subject only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generic filename / parameter details:
session_name{1} = 'CanonicalRS/';
nscans{1} = 55;
Fs_to_run{1} = 250;
nSj{1} = 55;
rawdatapath{1} = '/Volumes/CamsHD2/NottinghamRS/raw_data';

session_name{2} = 'Study1/';
nscans{2} = 42;
Fs_to_run{2} = [250,600];
nSj{2} = 21;
rawdatapath{2} = '/Volumes/CamsHD2/YunzheData/Replaydata4Cam/RawData';

session_name{3} ='Study2/';
nscans{3} = 44;
Fs_to_run{3} = [250,600];
nSj{3} = 22;
rawdatapath{3} = '/Volumes/CamsHD2/YunzheData/StrLearn_MEGexp/MEGData';

session_name{4} = 'Study1_FLI/';
nscans{4} = 2*21;
Fs_to_run{4} = 250;
nSj{4} = 21;

session_name{5} = 'Study2_FLI/';
nscans{5} = 3*22;
Fs_to_run{5} = 250;
nSj{5} = 22;

bestmodel=5;
K=12;
whichstudy=1

if whichstudy==1
  CanonicalRSN = false;
else
  CanonicalRSN = true; % canonical refers to RSNs trained on nottingham data and refitted
end

color_scheme = set1_cols();

% Define sample rate
sample_rate = 250;

if CanonicalRSN
  template_string = [int2str(bestmodel),'usingtemplate'];
else
  template_string = [int2str(bestmodel),''];
end

studydir = [wd,session_name{whichstudy},'250Hz/'];

savebase = fullfile( [wd,session_name{whichstudy},'Fig2_',template_string,'/' ])
if ~exist( savebase )
  mkdir(savebase);
end


hmmdir = [studydir,'hmm_1to45hz/'];
hmmfile = [hmmdir,'hmm',template_string,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm.mat'];
load(hmmfile);    % no need to permute here as will do below
prepdatafile = [hmmdir,'hmm_parc_giles_symmetric__pcdim80_voxelwise_embed14.mat'];
load(prepdatafile,'hmmT','subj_inds');
scan_T = cell2mat(hmmT);

parc_file = ['fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz'];
parc = parcellation(parc_file);
for k=1:hmm.K;statelabels{k}={'RSN-State ',int2str(k)};end

%% SECTION 1: Analyse Transition matrix and order states:

% Plot MDS map of state network:
disttoplot = plotMDS_states(hmm);
[~,new_state_ordering] = sort(disttoplot(:,1));

if any(new_state_ordering(2:end) ~= new_state_ordering(1:end-1)+1)
  hmm = hmm_permutestates(hmm,new_state_ordering);
  hmmfile = [hmmdir,'hmm',template_string,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm.mat'];
  save(hmmfile,'new_state_ordering','-append');
  disttoplot = disttoplot(new_state_ordering,:);
end
Gamma = hmm.gamma;

figure('Position', [440 519 391 279]);
for ik=1:hmm.K
  for k2=1:hmm.K
    if ik~=k2
      line([disttoplot(ik,1),disttoplot(k2,1)],[disttoplot(ik,2),disttoplot(k2,2)],...
        'color',0.5*[1,1,1]);hold on;
    end
  end
end
for ik=1:hmm.K
  scatter1 = scatter(disttoplot(ik,1),disttoplot(ik,2),400,...
    'MarkerFaceColor',color_scheme{ik},'MarkerEdgeColor',color_scheme{ik});
  hold on
  % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
  scatter1.MarkerFaceAlpha = 1;%.75;
  if ik==10
    scatter1 = scatter(disttoplot(ik-1,1),disttoplot(ik-1,2),400,...
      'MarkerFaceColor',color_scheme{ik-1},'MarkerEdgeColor',color_scheme{ik-1});
    hold on
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatter1.MarkerFaceAlpha = 0.5;
    text(disttoplot(ik-1,1)-0.03,disttoplot(ik-1,2),int2str(ik-1),'FontSize',12,'FontWeight','bold');hold on;
  end
  if ik<10
    text(disttoplot(ik,1)-0.03,disttoplot(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
  else
    text(disttoplot(ik,1)-0.05,disttoplot(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
  end
end
axis square
axis off

%% Look at each half of interstate intervals:
clear FO

dosimtest = false;

for iSj=1:max(hmm.subj_inds)
  fprintf(['\n Subject: ',int2str(iSj)]);
  vpath = Gamma(hmm.subj_inds==iSj,:) == repmat(max(Gamma(hmm.subj_inds==iSj,:),[],2),1,12);
  if dosimtest
    % replace vpath with simulated
    Ttemp = cumsum([hmm.P],2);
    statepath(1) = find(vpath(1,:));
    vpath(2:end,:) = zeros(size(vpath,1)-1,size(vpath,2));
    for iT=1:length(vpath)-1
      nextstate = find(rand(1)<Ttemp(find(vpath(iT,:)),:),1);
      statepath(iT+1) = nextstate;
      vpath(iT+1,nextstate) = 1;
    end
  end
  for ik=1:K
    tempaway = [];
    tempto = [];
    ontimes = find(diff([0;vpath(:,ik);0])==1);
    offtimes = find(diff([0;vpath(:,ik);0])==-1)-1;
    for t=1:length(offtimes)-1
      t_interval = ontimes(t+1)-offtimes(t);
      tempaway = [tempaway;vpath(offtimes(t):offtimes(t)+floor(t_interval/2),:)];
      tempto = [tempto;vpath(ontimes(t+1)-floor(t_interval/2):ontimes(t+1),:)];
      tintervals{iSj,ik}(t) = t_interval;
    end
    FO(ik,:,1,iSj) = mean(tempaway);
    FO(ik,:,2,iSj) = mean(tempto);
  end
end

%%

% be clear on dimensions: FO(ik1,ik2,:,:) refers to the concentration of
% state ik2 in the periods away from / towards state ik1

FO_m = mean(FO,4);
figure();imagesc(FO_m(:,:,1)-FO_m(:,:,2))

% ttests:
for ik1=1:12
  for ik2=1:12
    [~,pvals(ik1,ik2)] = ttest(squeeze(FO(ik1,ik2,1,:)-FO(ik1,ik2,2,:)));
  end
end
M = squeeze(mean(FO(:,:,1,:)-FO(:,:,2,:),4));
figure();imagesc(pvals<(0.05/144))

% conclusion: results are slightly less significant, but more
% interpretable?



%% plot with arrows:

disttoplot2 = disttoplot; % hack around with display:
disttoplot2(3,2) = -disttoplot2(3,2);
disttoplot2(9,2) = -disttoplot2(9,2);


sigpoints = pvals<(0.05);
figure('Position', [440 519 391 279]);
for ik=1:hmm.K
  for k2=1:hmm.K
    if sigpoints(ik,k2)
      %             line([disttoplot(ik,1),disttoplot(jk,1)],[disttoplot(ik,2),disttoplot(jk,2)],...
      %                 'color',0.5*[1,1,1]);hold on;
      if M(ik,k2)>0 % arrow from k1 to k2:
        quiver([disttoplot2(ik,1)],[disttoplot2(ik,2),],...
          disttoplot2(k2,1)-[disttoplot2(ik,1)],disttoplot2(k2,2)-[disttoplot2(ik,2)]);hold on;
      else % arrow from k2 to k1:
        quiver(disttoplot2(k2,1),disttoplot2(k2,2),...
          [disttoplot2(ik,1)]-disttoplot2(k2,1),[disttoplot2(ik,2)]-disttoplot2(k2,2));hold on;
      end
    end
  end
end
for ik=1:hmm.K
  scatter1 = scatter(disttoplot2(ik,1),disttoplot2(ik,2),400,...
    'MarkerFaceColor',color_scheme{ik},'MarkerEdgeColor',color_scheme{ik});
  hold on
  % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
  scatter1.MarkerFaceAlpha = 1;%.75;
  if ik==10
    scatter1 = scatter(disttoplot2(ik-1,1),disttoplot2(ik-1,2),400,...
      'MarkerFaceColor',color_scheme{ik-1},'MarkerEdgeColor',color_scheme{ik-1});
    hold on
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatter1.MarkerFaceAlpha = 0.5;
    text(disttoplot2(ik-1,1)-0.03,disttoplot2(ik-1,2),int2str(ik-1),'FontSize',12,'FontWeight','bold');hold on;
  end
  if ik<10
    text(disttoplot2(ik,1)-0.03,disttoplot2(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
  else
    text(disttoplot2(ik,1)-0.05,disttoplot2(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
  end
end
axis square
axis off

%% Determine geometric plot of network sequenceness for circular plot:
if 1
  % weights:
  beta = FO_m(:,:,1)-FO_m(:,:,2);
  
  % run permutations - note we can always specify 1 at same spot, and only
  % need consider 5 options for state 2 as others will be same but inverted
  myperms = perms(1:10);
  sequencemetric = zeros(length(myperms),5);
  for i2=1:length(myperms)
    if mod(i2,10000)==0
      fprintf(['\n Now up to run ',int2str(i2),' of ',int2str(size(myperms,1))]);
    end
    for i=1:5
      % setup state points on unit circle:
      %manualorder = [1,2,3,6,4,9,11,8,12,10,7,5];
      manualorder = [1,2+myperms(i2,:)];
      manualorder = [manualorder(1:i),2,manualorder(i+1:end)];
      disttoplot_manual = zeros(12,1);
      for i3=1:12
        disttoplot_manual(manualorder(i3)) = exp(sqrt(-1)*i3/12*2*pi);
      end
      angleplot = exp(sqrt(-1)*(angle(disttoplot_manual.')-angle(disttoplot_manual)));
      
      sequencemetric(i2,i) = imag(sum(sum(angleplot.*beta)));
    end
  end
  [~,m1] = max(max(abs(sequencemetric)));
  [trueseqmetric,m] = max(sequencemetric(:,m1));
  
  bestseq = [1,2+myperms(m,:)];
  bestseq = [bestseq(1:m1),2,bestseq(m1+1:end)];
  
  %bestseq = [1     2     3     4     6     8    11     9    12    10     7     5];
end

% we can also test whether the above 'true' sequenceness could arise by
% chance configuration of states - shuffle pairs of connection weights,
% determine the most sequential configuration that emerges and compare that
% to the non-shuffled value.

% this confirms this network pattern arises from the specific configuration
% of the entire network, not by chance plotting optimisation

for iperm=1:100
  fprintf(['\n Runnning perm ',int2str(iperm)]);
  % shuffle beta weights:
  %beta_perm = beta(randperm(12),:); % NO- should shuffle in pairs only
  beta_top = beta(triu(true(12),1));
  beta_low =beta';
  beta_low = beta_low(triu(true(12),1));
  indperm = randperm(66);
  inds = reshape(1:144,[12,12]);
  inds_top = inds(triu(true(12),1));
  inds_low = inds';
  inds_low = inds_low(triu(true(12),1));
  beta_perm = zeros(12);
  for i=1:length(beta_top)
    beta_perm(inds_top(i)) = beta_top(indperm(i));
    beta_perm(inds_low(i)) = beta_low(indperm(i));
  end
  
  % run permutations - note we can always specify 1 at same spot, and only
  % need consider 5 options for state 2 as others will be same but inverted
  myperms = perms(1:10);
  sequencemetric = zeros(length(myperms),5);
  for i2=1:length(myperms)
    %         if mod(i2,10000)==0
    %             fprintf(['\n Now up to run ',int2str(i2),' of ',int2str(size(myperms,1))]);
    %         end
    for i=1:5
      % setup state points on unit circle:
      %manualorder = [1,2,3,6,4,9,11,8,12,10,7,5];
      manualorder = [1,2+myperms(i2,:)];
      manualorder = [manualorder(1:i),2,manualorder(i+1:end)];
      disttoplot_manual = zeros(12,1);
      for i3=1:12
        disttoplot_manual(manualorder(i3)) = exp(sqrt(-1)*i3/12*2*pi);
      end
      angleplot = exp(sqrt(-1)*(angle(disttoplot_manual.')-angle(disttoplot_manual)));
      
      sequencemetric(i2,i) = imag(sum(sum(angleplot.*beta_perm)));
    end
  end
  [~,m_perm] = max(max(abs(sequencemetric)));
  [permseqmetric(iperm),~] = max(sequencemetric(:,m_perm));
  
end


%% now repeat as circular diagram with arrows:

figdir = ['/Users/chiggins/Google Drive/Doctoral Files/3.0 Experimental Work/6.0 RSN analysis/Study',int2str(whichstudy)];

%manualorder = [1,2,3,6,4,9,11,8,12,10,7,5];
manualorder = [1     2     3     4     6     8    11     9    12    10     7     5]

disttoplot_manual = zeros(12,2);
for i=1:12
  temp = exp(sqrt(-1)*(i+2)/12*2*pi);
  disttoplot_manual(manualorder(i),:) = [real(temp),imag(temp)];
end

sigpoints = pvals<(0.05/144);
figure('Position',[440 501 402 297]);
for ik1=1:hmm.K
  for k2=1:hmm.K
    if sigpoints(ik1,k2)
      %             line([disttoplot(ik,1),disttoplot(jk,1)],[disttoplot(ik,2),disttoplot(jk,2)],...
      %                 'color',0.5*[1,1,1]);hold on;
      linescale = sqrt(sum((disttoplot_manual(k2,:)-disttoplot_manual(ik1,:)).^2));
      if linescale>1.42 % ie line is four or more steps
        linescale = 1;
      elseif linescale >1.2% line is three steps
        linescale = 0.98;
      elseif linescale > 0.75 % line is two steps
        linescale = 0.9;
      else  % line is one step
        linescale = 0.8;
      end
      quivlength = linescale * sqrt(sum((disttoplot_manual(k2,:)-disttoplot_manual(ik1,:)).^2));
      if M(ik1,k2)>0 % arrow from k1 to k2:
        quiver([disttoplot_manual(ik1,1)],[disttoplot_manual(ik1,2),],...
          linescale*[disttoplot_manual(k2,1)-disttoplot_manual(ik1,1)],linescale*[disttoplot_manual(k2,2)-disttoplot_manual(ik1,2)],'Color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.4/quivlength);hold on;
      else % arrow from k2 to k1:
        quiver(disttoplot_manual(k2,1),disttoplot_manual(k2,2),...
          linescale*[disttoplot_manual(ik1,1)-disttoplot_manual(k2,1)],linescale*[disttoplot_manual(ik1,2)-disttoplot_manual(k2,2)],'Color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.4/quivlength);hold on;
      end
    end
  end
end
for ik=1:hmm.K
  scatter1 = scatter(disttoplot_manual(ik,1),disttoplot_manual(ik,2),400,...
    'MarkerFaceColor',color_scheme{ik},'MarkerEdgeColor',color_scheme{ik});
  hold on
  % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
  scatter1.MarkerFaceAlpha = 1;%.75;
  if ik==10
    scatter1 = scatter(disttoplot_manual(ik-1,1),disttoplot_manual(ik-1,2),400,...
      'MarkerFaceColor',color_scheme{ik-1},'MarkerEdgeColor',color_scheme{ik-1});
    hold on
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatter1.MarkerFaceAlpha = 0.5;
    text(disttoplot_manual(ik-1,1)-0.03,disttoplot_manual(ik-1,2),int2str(ik-1),'FontSize',12,'FontWeight','bold');hold on;
  end
  if ik<10
    text(disttoplot_manual(ik,1)-0.03,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
  else
    text(disttoplot_manual(ik,1)-0.05,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
  end
end
axis square
axis off
%title(['State ',int2str(ik1),' sequential patterns'])

print([figdir,'NetworkPlot'],'-dpng');

%% weight by individual state sequenceness:
beta = (FO_m(:,:,1)-FO_m(:,:,2))./mean(FO_m,3);
% compute state sequenceness metric:
disttoplot_manual = zeros(12,1);
for i3=1:12
  disttoplot_manual(manualorder(i3)) = exp(sqrt(-1)*i3/12*2*pi);
end
angleplot = exp(sqrt(-1)*(angle(disttoplot_manual.')-angle(disttoplot_manual)));

sequencemetricmat = imag(angleplot.*beta);
for i=1:12
  weights(i) = sum([sequencemetricmat(i,:),sequencemetricmat(:,i)']);
  %weights(i) = sum([sequencemetricmat(i,:)]);
  %weights(i) = sum([sequencemetricmat(:,i)']);
  %weights(i) = std(beta(i,:));
end
weights = weights./mean(weights);
manualorder = [1     2     3     4     6     8    11     9    12    10     7     5]
manualorder = [1,fliplr(manualorder(2:end))];
disttoplot_manual = zeros(12,2);
for i=1:12
  temp = exp(sqrt(-1)*(i+2)/12*2*pi);
  disttoplot_manual(manualorder(i),:) = weights(i)*[real(temp),imag(temp)];
end

sigpoints = pvals<(0.05/144);
figure('Position',[440 501 402 297]);
for ik1=1:hmm.K
  for k2=1:hmm.K
    if sigpoints(ik1,k2)
      %             line([disttoplot(ik,1),disttoplot(jk,1)],[disttoplot(ik,2),disttoplot(jk,2)],...
      %                 'color',0.5*[1,1,1]);hold on;
      linescale = sqrt(sum((disttoplot_manual(k2,:)-disttoplot_manual(ik1,:)).^2));
      if linescale>1.42 % ie line is four or more steps
        linescale = 1;
      elseif linescale >1.2% line is three steps
        linescale = 0.98;
      elseif linescale > 0.75 % line is two steps
        linescale = 0.9;
      else  % line is one step
        linescale = 0.8;
      end
      quivlength = linescale * sqrt(sum((disttoplot_manual(k2,:)-disttoplot_manual(ik1,:)).^2));
      if M(ik1,k2)>0 % arrow from k1 to k2:
        quiver([disttoplot_manual(ik1,1)],[disttoplot_manual(ik1,2),],...
          linescale*[disttoplot_manual(k2,1)-disttoplot_manual(ik1,1)],linescale*[disttoplot_manual(k2,2)-disttoplot_manual(ik1,2)],'Color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.4/quivlength);hold on;
      else % arrow from k2 to k1:
        quiver(disttoplot_manual(k2,1),disttoplot_manual(k2,2),...
          linescale*[disttoplot_manual(ik1,1)-disttoplot_manual(k2,1)],linescale*[disttoplot_manual(ik1,2)-disttoplot_manual(k2,2)],'Color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.4/quivlength);hold on;
      end
    end
  end
end
for ik=1:hmm.K
  scatter1 = scatter(disttoplot_manual(ik,1),disttoplot_manual(ik,2),400,...
    'MarkerFaceColor',color_scheme{ik},'MarkerEdgeColor',color_scheme{ik});
  hold on
  % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
  scatter1.MarkerFaceAlpha = 1;%.75;
  if ik==10
    scatter1 = scatter(disttoplot_manual(ik-1,1),disttoplot_manual(ik-1,2),400,...
      'MarkerFaceColor',color_scheme{ik-1},'MarkerEdgeColor',color_scheme{ik-1});
    hold on
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatter1.MarkerFaceAlpha = 0.5;
    text(disttoplot_manual(ik-1,1)-0.03,disttoplot_manual(ik-1,2),int2str(ik-1),'FontSize',12,'FontWeight','bold');hold on;
  end
  if ik<10
    text(disttoplot_manual(ik,1)-0.03,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
  else
    text(disttoplot_manual(ik,1)-0.05,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
  end
end
axis square
%axis off
% plot unit circle:
for t=1:1000
  unitcircle(t) = cos(2*pi*t/1000) + sqrt(-1)*sin(2*pi*t/1000);
end
plot(real(unitcircle),imag(unitcircle),'k--','LineWidth',1.5)

%title(['State ',int2str(ik1),' sequential patterns'])

print([figdir,'NetworkPlot_weighted'],'-dpng');


%% do figures per state:

disttoplot_manual = zeros(12,2);
for i=1:12
  temp = exp(sqrt(-1)*(i+2)/12*2*pi);
  disttoplot_manual(manualorder(i),:) = [real(temp),imag(temp)];
end

sigpoints = pvals<(0.05/144);
figure('Position',  [1 57 1440 741]);
for k2=1:hmm.K
  subplot(2,6,k2)
  
  for ik1=1:hmm.K
    if sigpoints(ik1,k2)
      %             line([disttoplot(ik,1),disttoplot(jk,1)],[disttoplot(ik,2),disttoplot(jk,2)],...
      %                 'color',0.5*[1,1,1]);hold on;
      linescale = sqrt(sum((disttoplot_manual(k2,:)-disttoplot_manual(ik1,:)).^2));
      if linescale>1.42 % ie line is four or more steps
        linescale = 1;
      elseif linescale >1.2% line is three steps
        linescale = 0.98;
      elseif linescale > 0.75 % line is two steps
        linescale = 0.9;
      else  % line is one step
        linescale = 0.8;
      end
      quivlength = linescale * sqrt(sum((disttoplot_manual(k2,:)-disttoplot_manual(ik1,:)).^2));
      if M(ik1,k2)>0 % arrow from k1 to k2:
        quiver([disttoplot_manual(ik1,1)],[disttoplot_manual(ik1,2),],...
          linescale*[disttoplot_manual(k2,1)-disttoplot_manual(ik1,1)],linescale*[disttoplot_manual(k2,2)-disttoplot_manual(ik1,2)],'Color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.4/quivlength);hold on;
      else % arrow from k2 to k1:
        quiver(disttoplot_manual(k2,1),disttoplot_manual(k2,2),...
          linescale*[disttoplot_manual(ik1,1)-disttoplot_manual(k2,1)],linescale*[disttoplot_manual(ik1,2)-disttoplot_manual(k2,2)],'Color',[0 0 0]+0.8,'LineWidth',2,'MaxHeadSize',0.4/quivlength);hold on;
      end
    end
  end
  
  for ik=1:hmm.K
    scatter1 = scatter(disttoplot_manual(ik,1),disttoplot_manual(ik,2),400,...
      'MarkerFaceColor',color_scheme{ik},'MarkerEdgeColor',color_scheme{ik});
    hold on
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatter1.MarkerFaceAlpha = 1;%.75;
    if ik==10
      scatter1 = scatter(disttoplot_manual(ik-1,1),disttoplot_manual(ik-1,2),400,...
        'MarkerFaceColor',color_scheme{ik-1},'MarkerEdgeColor',color_scheme{ik-1});
      hold on
      % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
      scatter1.MarkerFaceAlpha = 0.5;
      text(disttoplot_manual(ik-1,1)-0.03,disttoplot_manual(ik-1,2),int2str(ik-1),'FontSize',12,'FontWeight','bold');hold on;
    end
    if ik<10
      text(disttoplot_manual(ik,1)-0.03,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
    else
      text(disttoplot_manual(ik,1)-0.05,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
    end
  end
  axis square
  axis off
  t = title(['State ',int2str(k2),' sequential patterns'])
  p = t.Position .* [1,1.1,1];
  set(t,'Position',t.Position .* [1,1.2,1])
end
%print([figdir,'StateLeadLagAnalysis'],'-dpng');
print([figdir,'StateLeadLagAnalysis'],'-dsvg');

figure('Position', [440 579 114 219]);
quiver(0,0,1,0,'Color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.4/quivlength);hold on;
quiver(0,1,1,0,'Color',[0 0 0]+0.8,'LineWidth',2,'MaxHeadSize',0.4/quivlength);hold on;
axis off;
print([figdir,'StateLeadLagAnalysis_legend'],'-dpng');

%% toy plot to explain:

iSj=5%:max(hmm.subj_inds)
fprintf(['\n Subject: ',int2str(iSj)]);
vpath = Gamma(hmm.subj_inds==iSj,:);
[~,vpath] = max(vpath,[],2);
t_segment = 34300 + [1:720];

figure('Color','White','Position',[440 553 537 245]);
for k=1:12
  myline = NaN(length(t_segment),1);
  myline(vpath(t_segment)==k)=(13-k);
  plot((1/250):(1/250):(length(t_segment)/250),myline,'LineWidth',20,'Color',color_scheme{k});
  hold on;
  yaxislabels{13-k} = ['State ',int2str(k)];
end
%axis off;
set(gca,'YTick',[1:12]);
set(gca,'YTickLabels',yaxislabels);
grid off
h=gca;
set(gca,'xcolor','none')
h.XAxis.Label.Color=[0 0 0];
h.XAxis.Label.Visible='on';
xlabel('Time (sec)');
%set(gca,'ycolor','none')
h.YAxis.Label.Color=[0 0 0];
%h.YAxis.Label.Visible='on';
ylim([0.5,12.5])
print([figdir,'ExampleFigure1'],'-dpng');

%% Get coherence maps

for kk = 1:K
  figure();
  graph = abs(squeeze(nnmfWB_res.nnmf_coh_maps(kk,1,:,:)));
  tmp=squash(triu(graph));
  inds2=find(tmp>1e-10);
  data=tmp(inds2);
  
  S2=[];
  S2.data=squash(data);
  S2.do_fischer_xform=0;
  S2.pvalue_th=0.05/(38.^2);
  S2.do_plots=0;
  graph_ggm=teh_graph_gmm_fit(S2);
  
  th=graph_ggm.normalised_th;
  graph=graph_ggm.data';
  
  if th<1.96 % less than 2 stds from mean
    graph(graph<th)=NaN;
    graphmat=nan(nparcels, nparcels);
    graphmat(inds2)=graph;
    graph=graphmat;
  else
    % few sparse connections, do not plot:
    graph = nan(nparcels);
  end
  
  if all(isnan(graph(:)))
    graph(1,1)=1;
  end
  parcelFile = fullfile(osldir,'parcellations',parc_file);
  spatialRes = 8;
  spatialMap = nii.quickread(parcelFile, spatialRes);
  mni_coords = find_ROI_centres(spatialMap, spatialRes, 0, osldir);
  
  % plot
  nROIs = size(graph,2);
  colorLims = [th th+1];
  sphereCols = repmat([30 144 255]/255, nROIs, 1);
  edgeLims = [4 8];
  
  osl_braingraph(graph, colorLims, repmat(0.5,nROIs,1), [0 1], mni_coords, [], 0, sphereCols, edgeLims);
  
  view([0 90]);
  zoom(1);
  print([savebase '/WB_Coherence',int2str(kk),'_power'],'-depsc');
  
  viewZ = {[270,0],[-270,0],[0,90]};
  figure('Color', 'w','Position',[547 100 577 453]);
  ax(1) = axes('Position',[0 0.5 0.5 0.5]);
  ax(2) = axes('Position',[0.55 0.5 0.5 0.5]);
  ax(3) = axes('Position',[0.27 0.1 0.5 0.5]);
  
  % and plot 3 way brain graphs:
  for iplot=1:3
    axes(ax(iplot))
    osl_braingraph(graph, colorLims, repmat(0.5,nROIs,1), [0 1], mni_coords, [], 0, sphereCols, edgeLims);
    view(ax(iplot),viewZ{iplot})
    colorbar('hide')
  end
  print([savebase '/WB_Coherence',int2str(kk),'_3way'],'-depsc');
  close all;
end


%% make video of psd and coherence:

% plot metastate psd and coh:
studydir='/ohba/pi/mwoolrich/datasets/ReplayData/Neuron2020Analysis/CanonicalRS/250Hz/';
bestmodel=5;
template_string = [int2str(bestmodel),''];
K=12;
mtfilename  = [studydir, 'hmm_1to45hz/', 'hmm5_parc_giles_symmetric__pcdim80_voxelwise_embed14_K12_big1_dyn_modelhmm_store/state_netmats_mtsess_2_vn0_soft_global0.mat'];
hmmdir = [studydir, 'hmm_1to45hz/', 'hmm5_parc_giles_symmetric__pcdim80_voxelwise_embed14_K12_big1_dyn_modelhmm.mat'];
if ~exist(mtfilename)
  error('Soft state timecourses not found - rerun this analysis!');
end

[psd,coh] = loadMTspect(studydir,K,template_string);

% figdir = ['/Users/chiggins/Google Drive/Doctoral Files/3.0 Experimental Work/6.0 RSN analysis/Study',int2str(whichstudy)];
figdir=pwd;
%manualorder = [1,2,3,6,4,9,11,8,12,10,7,5];
usemanual=false;
manualorder = [1     2     3     4     6     8    11     9    12    10     7     5];
if usemanual
  bestseq = manualorder;
else
  optimalseqfile = [fileparts(hmmdir), '/bestseq',int2str(1),'.mat'];
  tmp=load(optimalseqfile);
  bestseq = tmp.bestsequencemetrics{1};  
end


disttoplot_manual = zeros(12,2);
for i=1:12
  temp = exp(sqrt(-1)*(i+2)/12*2*pi);
  disttoplot_manual(bestseq(i),:) = [real(temp),imag(temp)];
end


load([fileparts(hmmdir), '/HMMsummarymetrics.mat'])
FO = hmm_1stlevel.FO_intervals;
mean_direction = squeeze(mean(FO(:,:,1,:)-FO(:,:,2,:),4));
mean_assym = squeeze(mean((FO(:,:,1,:)-FO(:,:,2,:))./mean(FO,3),4));
for ik1=1:K
  for ik2=1:K
    pvals(ik1,ik2) = ttest(squeeze(FO(ik1,ik2,1,:)-FO(ik1,ik2,2,:)));
  end
end
sigpoints = pvals<(0.05/132);

hmm_1stlevel.FO_intervals = FO;
hmm_1stlevel.FO_assym = squeeze((FO(:,:,1,:)-FO(:,:,2,:))./mean(FO,3));


load(hmmdir)

nnmf_outfileWB=fullfile([studydir, 'hmm_1to45hz/','embedded_HMM_K',int2str(K),template_string,'_nnmfWB']);
load(nnmf_outfileWB)
nparcels=38;
Gamma = hmm.gamma;
net_mean = zeros(nparcels,size(Gamma,2));
thresh_mean = zeros(size(Gamma,2),1);
parc=parcellation('fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
for k = 1:size(Gamma,2)
  net_mean(:,k) = squeeze(nnmfWB_res.nnmf_psd_maps(k,1,:))';
end
for kk=1:hmm.K
  toplot = net_mean(:,kk);%-mean(net_mean,2);
  psdthresh = prctile(abs(toplot),50);
  %CL = psdthresh*[-1.25,1.25];
  CL = max(abs(squash(net_mean(:,:)))) * [0 1];
  %CL = [min(toplot(:)), max(toplot(:))];
  toplot(abs(toplot)<psdthresh) = NaN;
  %f2 = parc.plot_surface(psds,0,false,'enclosing',CL);
  f2 = plot_surface_4way(parc,toplot,0,false,'enclosing',[],[],CL);
  %    print([savebase '/St',int2str(kk),'_power_WB_4wayplot'],'-depsc');
  %    close(f2);
  colormap(cmap)
end

%% replay stuff
%
bestmodel=5;
    template_string = [int2str(bestmodel),'usingtemplate'];
 studydir='/ohba/pi/mwoolrich/datasets/ReplayData/Neuron2020Analysis/Study2/hmm_1to45hz/';
 hmmfile = [studydir,'hmm',template_string,'_parc_giles_symmetric__pcdim80_voxelwise_embed13_K',int2str(K),'_big1_dyn_modelhmm.mat'];
 load(hmmfile)
 disttoplot = plotMDS_states(hmm);
[~,new_state_ordering] = sort(disttoplot(:,1));

if any(new_state_ordering(2:end) ~= new_state_ordering(1:end-1)+1)
    hmm = hmm_permutestates(hmm,new_state_ordering);
    disttoplot = disttoplot(new_state_ordering,:);
end
 Gamma = hmm.gamma;
 
 datadir = '/ohba/pi/mwoolrich/datasets/ReplayData/Neuron2020Analysis/Study2/bfnew_1to45hz/'
 [maskA,maskB,triggerpoints,goodsamples] = getSubjectMasks(datadir);
 prepdatafile = [studydir,'hmm_parc_giles_symmetric__pcdim80_voxelwise_embed13.mat'];
load(prepdatafile,'hmmT','subj_inds');
scan_T = cell2mat(hmmT);
%  studydir='/ohba/pi/mwoolrich/datasets/ReplayData/Neuron2020Analysis/CanonicalRS/250Hz/';
% studydir = [wd,session_name{whichstudy},'250Hz/'];



load('/ohba/pi/mwoolrich/mvanes/Projects/Replay/Study2/STUDYII_PreplayOnset.mat')
replayScores(1:2:(2*22),:) = ToRall;
load('/ohba/pi/mwoolrich/mvanes/Projects/Replay/Study2/STUDYII_ReplayOnset.mat')
replayScores(2:2:(2*22),:) = ToRall;

t_window=125;
R = [[1;1+cumsum(scan_T(1:end-1)')],[cumsum(scan_T(1:end)')]];
if ~all(R(2:end,1)==R(1:end-1,2)+1)
    % correct any uncorrected offsets in R:
    R(2:end,1) = R(1:end-1,2)+1;
end
for iSj=1:44
    Gam_long{iSj} = NaN(length(goodsamples{iSj}),K);
    Gam_long{iSj}(goodsamples{iSj},:) = Gamma(R(iSj,1):R(iSj,2),:);
    Gam_long{iSj} = Gam_long{iSj}(triggerpoints{iSj},:);
    replaytimes = replayScores(iSj,:) > prctile(replayScores(iSj,:),99);
    %eliminate adjacent points
    replaytimes= [0,diff(replaytimes)]==1;
    %eliminate border points:
    replaytimes(1:t_window)=0;replaytimes([end-t_window]:end)=0;
    t_i{iSj} = round(find(replaytimes)*2.5);
    [~, vpath{iSj}] = max(Gam_long{iSj},[],2, 'omitnan');
    
    %         if length(Gam_long)>length(replaytimes)
    %         Q = length(Gam_long) / length(replaytimes);
    %         t_g{iSj} = round(t_i{iSj}*Q);
    %     else
    %         Q = 1;
    %         t_g{iSj} = t_i;
    %     end
%     
%         if any(t_g{iSj}>size(Gam_long,1)-t_window)
%         t_g{iSj}(t_g{iSj}>size(Gam_long,1)-t_window) = [];
%         t_i{iSj}(t_i{iSj}>size(Gam_long,1)-t_window) = [];
%         end
%     
%             Gamma_mu = nanmean(Gam_long);
%     X_subj = Gam_long-repmat(Gamma_mu,length(Gam_long),1);
%     Gammat = NaN(2*t_window+1,K,length(t_i));
%     for iEv=1:length(t_i)
%         t_sample_g = [t_g(iEv)-t_window : t_g(iEv)+t_window];
%         Gammat(:,:,iEv) = X_subj(t_sample_g,:);
%         
%     end
end

% actualclock = [1,2,5,7,10,12,9,11,8,6,4,3];

for k=1:44
%   v2{k} = [];
  tmp = vpath{k}(t_i{k});
  for ik=1:12
    replay_counts(k,ik) = sum(tmp==ik);
%     v2{k} = [v2{k}; disttoplot_manual(ik,1)*ones(replay_counts(k,ik),1)];
  end
end

for k=1:44, Q{k} = 75000; end
[FO,pvals,t_intervals] = computeLongTermAsymmetry(vpath, Q,K);
bonf_ncomparisons = K.^2-K;
mean_direction = squeeze(mean(FO(:,:,1,:)-FO(:,:,2,:),4));
mean_assym = squeeze(mean((FO(:,:,1,:)-FO(:,:,2,:))./mean(FO,3),4));
sigpoints = pvals<(0.05);
%   
%   v{k}=zeros(75000,1);
%   for i=1:12
%     v{k}(vpath{k}==i) = disttoplot_manual(i,1)+sqrt(-1)*disttoplot_manual(i,2);
%   end
% 
% for k=1:44
%   Q{k} = v{k}(t_i{k});
% end
% 
% percreplay=zeros(12,22);
% for k=1:22
%   for i=1:12
%     percreplay(i,k)=sum(Q{k}==v{k}(find(vpath{k}==i,1)))/length(Q{k});
%   end
% end
% percreplay_u = mean(percreplay,2);
% percreplay_sem = std(percrepfiglay,[],2)/sqrt(22);
% 
% 
% for k=1:44
%   Q2(k) = circ_mean(angle(v2{k}));
% end
%%
if doreplay
  V=VideoWriter('/ohba/pi/mwoolrich/mvanes/scripts/Tinda/TINDA_replay.avi');
else
  
end
V.FrameRate=7.5;
open(V)
for ivid=[numframes/4:-1:1, numframes:-1:numframes/4]
  %%%%%%%%%
  % CLOCK %
  %%%%%%%%%
fig=figure('Position',[0 0 570 413]);

  temp = exp(sqrt(-1)*(ivid)/numframes*2*pi);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CREATE MOVING AVERAGE WEIGHTS %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  q = abs((disttoplot_manual) - temp);
  [val, idx]=sort(abs(q), 'ascend');
  xpdf = zeros(12,1);
  if val(1)<0.001
    xpdf(idx(1))=1;
  else
    xpdf(idx(1:2)) = 1./val(1:2);
  end
  xpdf = xpdf./sum(xpdf);

  

  dum=raincloud_plot(v*xpdf, 'color', [0.8500 0.3250 0.0980]), xlim([0 60]); ylim([-0.1 0.12]),axis off,
  
  f=getframe(fig);
  size(f)
  writeVideo(V,f)
  close(fig)
end

close(V)
%% make TINDA movie
doreplay=1;
color_scheme = set1_cols();
cmap = colormap('inferno');
numframes=120;%480;
if isfile([figdir,'/tinda_full_4.avi'])
  delete([figdir,'/tinda_full_4.avi'])
end
if doreplay
  v=VideoWriter('/ohba/pi/mwoolrich/mvanes/scripts/Tinda/TINDA_replay.avi');
else
  v=VideoWriter('/ohba/pi/mwoolrich/mvanes/scripts/Tinda/TINDA.avi');
end
v.FrameRate = 7.5;
open(v)
nrow = 6;
ncol = 6;

% some coherence topo set up
parcelFile = fullfile(osldir,'parcellations',    'fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
spatialRes = 8;
spatialMap = nii.quickread(parcelFile, spatialRes);
mni_coords = find_ROI_centres(spatialMap, spatialRes, 0, osldir);
viewZ = {[270,0],[-270,0],[0,90]};

for ivid=[numframes/4:-1:1, numframes:-1:numframes/4]
  %%%%%%%%%
  % CLOCK %
  %%%%%%%%%
  fig=figure('Position',2*[0 0 560 420]);

  ax(5)=subplot(nrow,ncol,[3,4,9,10]);
%   for ik=1:hmm.K
%     scatter1 = scatter(disttoplot_manual(ik,1),disttoplot_manual(ik,2),400,...
%       'MarkerFaceColor',color_scheme{ik},'MarkerEdgeColor',color_scheme{ik});
%     hold on
%     % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
%     scatter1.MarkerFaceAlpha = 1;%.75;
%     if ik==10
%       scatter1 = scatter(disttoplot_manual(ik-1,1),disttoplot_manual(ik-1,2),400,...
%         'MarkerFaceColor',color_scheme{ik-1},'MarkerEdgeColor',color_scheme{ik-1});
%       hold on
%       % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
%       scatter1.MarkerFaceAlpha = 0.5;
%       text(disttoplot_manual(ik-1,1)-0.05,disttoplot_manual(ik-1,2),int2str(ik-1),'FontSize',12,'FontWeight','bold');hold on;
%     end
%     if ik<10
%       text(disttoplot_manual(ik,1)-0.05,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
%     else
%       text(disttoplot_manual(ik,1)-0.1,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
%     end
%   end
%   axis square
%   axis off
  cyclicalstateplot(bestseq,mean_direction,sigpoints,[],0)
  %plot vector simulated:
  hold on;
  temp = exp(sqrt(-1)*(ivid)/numframes*2*pi);
  plot([0 real(temp)],[0,imag(temp)],'k','LineWidth',2);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CREATE MOVING AVERAGE WEIGHTS %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %     xpdf = mvnpdf(disttoplot_manual,[real(temp),imag(temp)],0.5*eye(2));
%   q = abs((disttoplot_manual(:,1) + sqrt(-1)*disttoplot_manual(:,2)) - temp);
  q = abs((disttoplot_manual) - temp);
  [val, idx]=sort(abs(q), 'ascend');
  xpdf = zeros(12,1);
  if val(1)<0.001
    xpdf(idx(1))=1;
  else
    xpdf(idx(1:2)) = 1./val(1:2);
  end
  xpdf = xpdf./sum(xpdf);
  if doreplay
    %%%%%%%%%%%%%%%
    % REPLAY PLOT %
    %%%%%%%%%%%%%%%
    ax(6)=subplot(nrow,ncol, [5,11]);
    bar(sum(replay_counts*xpdf)), ylim([0 1200]), axis off, title('Replay counts'), 
    set(gca,'fontsize',16);
  end
  
  %%%%%%%%%%%%
  % PSD PLOT %
  %%%%%%%%%%%%
  I = logical(eye(38));
  psd_st = mean(psd(:,:,:,I),4);
  psd_temp = squeeze(sum(psd_st.*repmat(xpdf',[55,1,90]),2));
  subplot(nrow,ncol,[15,16,21,22])%11);
  shadedErrorBar(linspace(0,50,90),mean(psd_temp,1),std(psd_temp,[],1)./sqrt(55),{'-k', 'LineWidth', 2});
  q1=mean(psd_st) + std(psd_st);
  yl=max(q1(:));
  ylim([0,yl]);xlim([0,50]);
%   xlabel('Freq (Hz)')
  ylabel('PSD')

  %%%%%%%%%%%%%%%%%%
  % COHERENCE PLOT %
  %%%%%%%%%%%%%%%%%%
  coh_st = mean(coh(:,:,:,~I),4);
  coh_temp = squeeze(sum(coh_st.*repmat(xpdf',[55,1,90]),2));
  subplot(nrow,ncol,[27,28,33,34])%12);
  shadedErrorBar(linspace(0,50,90),mean(coh_temp,1),std(coh_temp,[],1)./sqrt(55),{'-k', 'LineWidth', 2});
  q1=mean(coh_st) + std(coh_st);
  yl=max(q1(:));
  ylim([0,yl]);xlim([0,50]);
  ylim([0.05,yl]);
  xlim([0,50]);
  xlabel('Freq (Hz)')
  ylabel('Coherence')
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  % COHERENCE TOPOGRAPHY %
  %%%%%%%%%%%%%%%%%%%%%%%%
  graph = abs(reshape(xpdf'*reshape(squeeze(nnmfWB_res.nnmf_coh_maps(:,1,:,:)), K, []), 38, 38));
  tmp=squash(triu(graph));
  inds2=find(tmp>1e-10);
  data=tmp(inds2);
  S2=[];
  S2.data=squash(data);
  S2.do_fischer_xform=0;
  S2.pvalue_th=0.05/(38.^2);
  S2.do_plots=0;
  graph_ggm=teh_graph_gmm_fit(S2);
  
  th=graph_ggm.normalised_th;
  graph=graph_ggm.data';
  if th<1.96 % less than 2 stds from mean
    graph(graph<th)=NaN;
    graphmat=nan(nparcels, nparcels);
    graphmat(inds2)=graph;
    graph=graphmat;
  else
    % few sparse connections, do not plot:
    graph = nan(nparcels);
  end
  if all(isnan(graph(:)))
    graph(1,1)=1;
  end
  
  nROIs = size(graph,2);
  colorLims = [th th+1];
  sphereCols = repmat([30 144 255]/255, nROIs, 1);
  edgeLims = [4 8];
  
%   figure('Color', 'w','Position',[547 100 577 453]);
  ax(3) = subplot(nrow,ncol,[25,26,31,32]);%axes('Position',[0 0.5 0.5 0.5]);
  ax(4) = subplot(nrow,ncol,[29,30,35,36]);%axes('Position',[0.55 0.5 0.5 0.5]);
%   ax(5) = %axes('Position',[0.27 0.1 0.5 0.5]);
  
  % and plot 3 way brain graphs:
  for iplot=1:2%3
    axes(ax(2+iplot))
    osl_braingraph(graph, colorLims, repmat(0.5,nROIs,1), [0 1], mni_coords, [], 0, sphereCols, edgeLims);
    view(ax(2+iplot),viewZ{iplot})
    colorbar('hide')
  end
    
  %%%%%%%%%%%%%%%%%%%%
  % POWER TOPOGRAPHY %
  %%%%%%%%%%%%%%%%%%%%
  toplot = net_mean*xpdf;%-mean(net_mean,2);
  if 1 % scale colorlim per power map
    CL = max(abs(squash(toplot(:,:)))) * [0 1];
  else % scale color lim over all power maps
    CL = max(abs(squash(net_mean(:,:)))) * [0 1];
  end
  if 0 % threshold the power map
    psdthresh = prctile(abs(toplot),50);
    toplot(abs(toplot)<psdthresh) = NaN;
  end

  ax(1)=subplot(nrow,ncol,[13,14,19,20]); axis off
  ax(2)=subplot(nrow, ncol, [17,18,23,24]); axis off

  plot_surface(parc,toplot,0,false,[], ax(1:2), CL);
  colormap(cmap)
  
  %%%%%%%%%%%%%%%
  % WRITE VIDEO %
  %%%%%%%%%%%%%%%
  f=getframe(fig);
  writeVideo(v,f)
  close(gcf)
end
close(v);

%% compare adjacent spatial correlation vs non adjacent:
temp = squeeze(mean(abs(psd),1));
I = logical(eye(38));
clear C
for ik1=1:12
  for ik2=1:12;
    for ifreq=1:90
      C(ik1,ik2,ifreq) = corr(squeeze(temp(ik1,ifreq,I)),squeeze(temp(ik2,ifreq,I)));
    end
  end
end
f = linspace(0,50,90);
figure();imagesc(C(:,:,find(f>10,1)));

%% check mean adjacent
adjacenymatrix = false(12);
adjacenymatrix(manualorder(1),manualorder(12)) = true;
for ik1=2:12;
  adjacenymatrix(manualorder(ik1),manualorder(ik1-1)) = true;
end
adjacenymatrix = logical(adjacenymatrix + adjacenymatrix');
C2 = permute(C,[3,1,2]);
for ifreq = 1:90
  cvals(ifreq,1) = mean(C2(ifreq,adjacenymatrix));
  cvals(ifreq,2) = mean(C2(ifreq,logical(~adjacenymatrix - eye(12))));
end
figure();plot(f,cvals,'LineWidth',2)
title('Adjacent states have higher spatial correlation than non-adjacent');
plot4paper('Freq (Hz)','Correlation');
legend('Adjacent','Non-adjacent');

%% redo for percentile length of intervals:


dosimtest = false;

runonquartiles = true;

intervalpercentiles = [0,20,40,60,80,100];
FO_prctile = zeros(K,K,2,max(hmm.subj_inds),length(intervalpercentiles)-1);
intervaltracker = cell(5,K);
for iSj=1:max(hmm.subj_inds)
  fprintf(['\n Subject: ',int2str(iSj)]);
  vpath = Gamma(hmm.subj_inds==iSj,:) == repmat(max(Gamma(hmm.subj_inds==iSj,:),[],2),1,12);
  if dosimtest
    % replace vpath with simulated
    Ttemp = cumsum([hmm.P],2);
    statepath(1) = find(vpath(1,:));
    vpath(2:end,:) = zeros(size(vpath,1)-1,size(vpath,2));
    for iT=1:length(vpath)-1
      nextstate = find(rand(1)<Ttemp(find(vpath(iT,:)),:),1);
      statepath(iT+1) = nextstate;
      vpath(iT+1,nextstate) = 1;
    end
  end
  for ik=1:K
    
    ontimes = find(diff([0;vpath(:,ik);0])==1);
    offtimes = find(diff([0;vpath(:,ik);0])==-1)-1;
    tintervals{iSj,ik} = ontimes(2:end)-offtimes(1:end-1);
    for ip=1:length(intervalpercentiles)-1
      p_low = prctile(tintervals{iSj,ik},intervalpercentiles(ip));
      p_high = prctile(tintervals{iSj,ik},intervalpercentiles(ip+1));
      tempaway = [];
      tempto = [];
      for t=1:length(offtimes)-1
        t_interval = ontimes(t+1)-offtimes(t);
        if t_interval>=p_low && t_interval<=p_high
          if ~runonquartiles
            tempaway = [tempaway;vpath(offtimes(t):offtimes(t)+floor(t_interval/2),:)];
            tempto = [tempto;vpath(ontimes(t+1)-floor(t_interval/2):ontimes(t+1),:)];
          else
            tempaway = [tempaway;vpath(offtimes(t):offtimes(t)+floor(t_interval/4),:)];
            tempto = [tempto;vpath(ontimes(t+1)-floor(t_interval/4):ontimes(t+1),:)];
          end
          intervaltracker{ip,ik}(length(intervaltracker{ip,ik})+1) = t_interval;
        end
      end
      FO_prctile(ik,:,1,iSj,ip) = mean(tempaway);
      FO_prctile(ik,:,2,iSj,ip) = mean(tempto);
    end
  end
end

% ttests:
pvals_prctile = zeros(K,K,length(intervalpercentiles)-1);
for ip=1:length(intervalpercentiles)-1
  for ik1=1:12
    for ik2=1:12
      [~,pvals_prctile(ik1,ik2,ip)] = ttest(squeeze(FO_prctile(ik1,ik2,1,:,ip)-FO_prctile(ik1,ik2,2,:,ip)));
    end
  end
end
M_prctile = squeeze(mean(FO_prctile(:,:,1,:,:)-FO_prctile(:,:,2,:,:),4));
figure();
for ip=1:length(intervalpercentiles)-1
  subplot(3,2,ip)
  imagesc(pvals_prctile(:,:,ip)<(0.05/144))
end


%% and plot for each percentile:

figdir = ['/Users/chiggins/Google Drive/Doctoral Files/3.0 Experimental Work/6.0 RSN analysis/'];

%manualorder = [1,2,3,6,4,9,11,8,12,10,7,5];
manualorder = [1     2     3     4     6     8    11     9    12    10     7     5]

disttoplot_manual = zeros(12,2);
for i=1:12
  temp = exp(sqrt(-1)*(i+2)/12*2*pi);
  disttoplot_manual(manualorder(i),:) = [real(temp),imag(temp)];
end
figure('Position', [1 56 1440 742]);
for ip=1:length(intervalpercentiles)-1
  sigpoints = pvals_prctile(:,:,ip)<(0.05/144);
  subplot(3,5,ip);
  for ik=1:hmm.K
    for k2=1:hmm.K
      
      if sigpoints(ik,k2)
        linescale = sqrt(sum((disttoplot_manual(k2,:)-disttoplot_manual(ik,:)).^2));
        if linescale>1.42 % ie line is four or more steps
          linescale = 1;
        elseif linescale >1.2% line is three steps
          linescale = 0.98;
        elseif linescale > 0.75 % line is two steps
          linescale = 0.9;
        else  % line is one step
          linescale = 0.8;
        end
        quivlength = linescale * sqrt(sum((disttoplot_manual(k2,:)-disttoplot_manual(ik,:)).^2));
        if M_prctile(ik,k2,ip)>0 % arrow from k1 to k2:
          quiver([disttoplot_manual(ik,1)],[disttoplot_manual(ik,2),],...
            linescale*[disttoplot_manual(k2,1)-disttoplot_manual(ik,1)],linescale*[disttoplot_manual(k2,2)-disttoplot_manual(ik,2)],'Color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.4/quivlength);hold on;
        else % arrow from k2 to k1:
          quiver(disttoplot_manual(k2,1),disttoplot_manual(k2,2),...
            linescale*[disttoplot_manual(ik,1)-disttoplot_manual(k2,1)],linescale*[disttoplot_manual(ik,2)-disttoplot_manual(k2,2)],'Color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.4/quivlength);hold on;
        end
        %             line([disttoplot(ik,1),disttoplot(jk,1)],[disttoplot(ik,2),disttoplot(jk,2)],...
        %                 'color',0.5*[1,1,1]);hold on;
        %                 if M_prctile(ik,k2,ip)>0 % arrow from k1 to k2:
        %                       quiver([disttoplot_manual(ik,1)],[disttoplot_manual(ik,2),],...
        %                           0.75*[disttoplot_manual(k2,1)-disttoplot_manual(ik,1)],0.85*[disttoplot_manual(k2,2)-disttoplot_manual(ik,2)],'Color',[0 0 0]);hold on;
        %                 else % arrow from k2 to k1:
        %                     quiver(disttoplot_manual(k2,1),disttoplot_manual(k2,2),...
        %                         0.75*[disttoplot_manual(ik,1)-disttoplot_manual(k2,1)],0.85*[disttoplot_manual(ik,2)-disttoplot_manual(k2,2)],'Color',[0 0 0]);hold on;
        %                 end
      end
    end
  end
  for ik=1:hmm.K
    scatter1 = scatter(disttoplot_manual(ik,1),disttoplot_manual(ik,2),400,...
      'MarkerFaceColor',color_scheme{ik},'MarkerEdgeColor',color_scheme{ik});
    hold on
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatter1.MarkerFaceAlpha = 1;%.75;
    if ik==10
      scatter1 = scatter(disttoplot_manual(ik-1,1),disttoplot_manual(ik-1,2),400,...
        'MarkerFaceColor',color_scheme{ik-1},'MarkerEdgeColor',color_scheme{ik-1});
      hold on
      % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
      scatter1.MarkerFaceAlpha = 0.5;
      text(disttoplot_manual(ik-1,1)-0.03,disttoplot_manual(ik-1,2),int2str(ik-1),'FontSize',12,'FontWeight','bold');hold on;
    end
    if ik<10
      text(disttoplot_manual(ik,1)-0.03,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
    else
      text(disttoplot_manual(ik,1)-0.05,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
    end
  end
  axis square
  axis off
  
  titlelabels = {'Shortest 20% intervals','20th to 40th percentile of intervals','40th to 60th percentile intervals','60th to 80th percentile of intervals','Longest 20% of intervals'};
  T = title(titlelabels{ip});
  T.Position(2)=1.5;
  %set(T,
  subplot(3,5,5+ip);
  distributionPlot(intervaltracker(ip,:),'showMM',2,'color',{color_scheme{1:size(FO,2)}})
  title('Interval Times');plot4paper('State','Time (ms)');grid on;
  %YL = 1.1*max(LTmerged(:))./ sample_rate * 1000;
  %set(gca,'YLim',[0 YL],'FontSize',fontsize,'FontSize',fontsize);
  % convert YTick to right scale:
  %m = cellfun(@prctile,intervaltracker(ip,:),90);
  m = prctile([intervaltracker{ip,:}],99.5);
  ylim([0,m]);
  YT = get(gca,'YTick');
  YT = YT./ sample_rate * 1000;
  set(gca,'YTickLabel',YT);
end
M = mean(FO_prctile,4);
M = squeeze(mean(M,3));
%figure('Position',[9 473 1432 325]);
% for ip=1:5
%     subplot(3,5,10+ip);
%     distancevector = M(:,:,ip)*disttoplot_manual;
%     distancevector = distancevector*1;
%     for ik=1:hmm.K
%         scatter1 = scatter(disttoplot_manual(ik,1),disttoplot_manual(ik,2),400,...
%             'MarkerFaceColor',color_scheme{ik},'MarkerEdgeColor',color_scheme{ik});
%         hold on
%         % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
%         scatter1.MarkerFaceAlpha = 1;%.75;
%         if ik==10
%             scatter1 = scatter(disttoplot_manual(ik-1,1),disttoplot_manual(ik-1,2),400,...
%             'MarkerFaceColor',color_scheme{ik-1},'MarkerEdgeColor',color_scheme{ik-1});
%             hold on
%             % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
%             scatter1.MarkerFaceAlpha = 0.5;
%             text(disttoplot_manual(ik-1,1)-0.03,disttoplot_manual(ik-1,2),int2str(ik-1),'FontSize',12,'FontWeight','bold');hold on;
%         end
%         if ik<10
%             text(disttoplot_manual(ik,1)-0.03,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
%         else
%             text(disttoplot_manual(ik,1)-0.05,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
%         end
%         quiver(0,0,distancevector(ik,1),distancevector(ik,2),'Color',color_scheme{ik},'LineWidth',2);hold on;
%
%     end
%     axis square
%     axis off
% end
if ~runonquartiles
  print([figdir,'cyclicpatterns_percentiled'],'-dpng')
else
  print([figdir,'cyclicpatterns_percentiled_runonquartiles'],'-dpng')
end

%% Plot the state density for each percentile of interval:

%FO_prctile(ik,:,1,iSj,ip)

M = mean(FO_prctile,4);
M = squeeze(mean(M,3));
figure('Position',[9 473 1432 325]);
for ip=1:5
  subplot(1,5,ip);
  distancevector = M(:,:,ip)*disttoplot_manual;
  distancevector = distancevector*1;
  for ik=1:hmm.K
    scatter1 = scatter(disttoplot_manual(ik,1),disttoplot_manual(ik,2),400,...
      'MarkerFaceColor',color_scheme{ik},'MarkerEdgeColor',color_scheme{ik});
    hold on
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatter1.MarkerFaceAlpha = 1;%.75;
    if ik==10
      scatter1 = scatter(disttoplot_manual(ik-1,1),disttoplot_manual(ik-1,2),400,...
        'MarkerFaceColor',color_scheme{ik-1},'MarkerEdgeColor',color_scheme{ik-1});
      hold on
      % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
      scatter1.MarkerFaceAlpha = 0.5;
      text(disttoplot_manual(ik-1,1)-0.03,disttoplot_manual(ik-1,2),int2str(ik-1),'FontSize',12,'FontWeight','bold');hold on;
    end
    if ik<10
      text(disttoplot_manual(ik,1)-0.03,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
    else
      text(disttoplot_manual(ik,1)-0.05,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
    end
    quiver(0,0,distancevector(ik,1),distancevector(ik,2),'Color',color_scheme{ik},'LineWidth',2);hold on;
    
  end
  axis square
  axis off
end
if ~runonquartiles
  print([figdir,'cyclicpatterns_percentiled_vector'],'-dpng')
else
  print([figdir,'cyclicpatterns_percentiled_runonquartiles_vector'],'-dpng')
end