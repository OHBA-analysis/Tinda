function [num_transitions,microsequences] = getMicroStateCycleVisits(vpath,T,vpathMicroStates)

if size(vpath,2)>1
    % convert from Gamma to vpath:
    [~,vpath] = max(vpath,[],2);
end

num_transitions = [];
for i=1:length(T)
    seq = vpath((sum(T(1:i-1))+1):sum(T(1:i)));
    seq_2ndlevel = vpathMicroStates((sum(T(1:i-1))+1):sum(T(1:i)));
    firstchange = find(diff(seq),1);
    num_trans = [];
    if ~isempty(firstchange)
        fullcycles = find(seq(1:end-1)==seq(firstchange) & seq(2:end)==seq(firstchange+1));
        %num_fullcycles = num_fullcycles + length(fullcycles);
        %time_tot = time_tot + fullcycles(end)-fullcycles(1);
        %cycletimes = [cycletimes;diff(fullcycles)];

        % fullcycles are the times at which a cycle is complete:
        for i2=1:length(fullcycles)-1
            microsequence = seq_2ndlevel(fullcycles(i2):fullcycles(i2+1));
            num_trans(i2,1) = sum(diff(microsequence)~=0); % total number of state transitions
            microsequence = seq_2ndlevel((fullcycles(i2)+1):fullcycles(i2+1));
            num_trans(i2,2) = length(unique(microsequence)); % total number of UNIQUE states observed
            microsequences{i2,i} = unique(microsequence);
        end
    end
    num_transitions = [num_transitions;num_trans];
end
if length(T)>1
    emptypoints = cellfun(@isempty,microsequences);
    microsequences = microsequences(~emptypoints);
end