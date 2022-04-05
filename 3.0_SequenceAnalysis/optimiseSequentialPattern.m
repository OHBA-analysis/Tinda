function bestseq = optimiseSequentialPattern(FO,optimalseqfile)
% this function reads in the mean pattern of differential fractional
% occupancy and computes the optimial display for a sequential circular
% plot visualisation.
%
% we test different possible metrics:
% sequence metric 1 is the group average direction. 
% sequence metric 2 is the average subject assymmetry (i.e. average of the
% normalized subject direction)
% sequence metric 3 is the group average assymmetry (i.e. metric 1
% normalized by the mean over all subjects

metric{1} = squeeze(mean(FO(:,:,1,:)-FO(:,:,2,:),4));
metric{2} = squeeze(mean((FO(:,:,1,:)-FO(:,:,2,:))./mean(FO,3),4));
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
        
        % metric 2 maximises the one-step (ie local) FO difference:
        for i3=1:3
            temp = 0;
            for k=1:11
                temp = temp + metric{i3}(manualorder(k),manualorder(k+1)) - metric{i3}(manualorder(k+1),manualorder(k));
            end
            sequencemetric{3+i3}(i2,i) = sum(temp);
        end
        
        % metric 3 
    end
end
for i=1:9
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