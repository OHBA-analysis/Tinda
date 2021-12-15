function bestseq = optimiseSequentialPattern(FO,optimalseqfile)
% this function reads in the mean pattern of differential fractional
% occupancy and computes the optimial display for a sequential circular
% plot visualisation.
%
% we test different possible metrics:
% sequence metric 1 is the mean FO assymetry
% sequence metric 2 is the proportional FO assymetry (ie assymetry as a
% proportion of a baseline - which is time spent in the state)
% sequence metric 3 is the propotional FO assymetry using the global,
% rather than subject-sepcific, baseline FO
%

metric{1} = squeeze(mean(FO(:,:,1,:)-FO(:,:,2,:),4));
temp = (FO(:,:,1,:)-FO(:,:,2,:))./mean(FO,3); % control for NANs where a state isn't visited
temp(mean(FO,3)==0)=0;
metric{2} = squeeze(mean(temp,4));
metric{3} = squeeze(mean(FO(:,:,1,:)-FO(:,:,2,:),4)) ./ mean(FO(:,:,:),3);

myperms = perms(1:10);
for i=1:9
    sequencemetric{i} = zeros(length(myperms),5);
end
for i2=1:length(myperms)
    if mod(i2,10000)==0
        fprintf(['\n Now up to run ',int2str(i2),' of ',int2str(size(myperms,1))]);
    end
    for i=1:5
        % setup state points on unit circle:
        manualorder = [1,2+myperms(i2,:)];
        manualorder = [manualorder(1:i),2,manualorder(i+1:end)];
        disttoplot_manual = zeros(12,1);
        for i3=1:12
            disttoplot_manual(manualorder(i3)) = exp(sqrt(-1)*i3/12*2*pi);
        end
        angleplot = exp(sqrt(-1)*(angle(disttoplot_manual.')-angle(disttoplot_manual)));
        
        for i3=1:3  
            sequencemetric{i3}(i2,i) = imag(sum(sum(angleplot.*metric{i3})));
        end
        
       
    end
end
for i=1:3
    [~,m1] = max(max(abs(sequencemetric{i})));
    [~,m] = max(abs(sequencemetric{i}(:,m1)));
    trueseqmetric = sequencemetric{i}(m,m1);

    bestseq{i} = [1,2+myperms(m,:)];
    bestseq{i} = [bestseq{i}(1:m1),2,bestseq{i}(m1+1:end)];
    % want rotation to be clockwise, so flip if necessary:
    if trueseqmetric>0
        bestseq{i} = [1,fliplr(bestseq{i}(2:end))];
    end
end
end