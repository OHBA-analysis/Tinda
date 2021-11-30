function vpath_new = simulateVpath(vpath,T,K)
% this function takes an inferred viterbi path and simulates another using
% only the transition probability structure expressed in the original
% viterbi path.

vpath_new = zeros(size(vpath));
if sum(T)~=length(vpath)
    error('T and vpath not same dimension');
end
% infer FO
for ik=1:K
    FO(ik) = mean(vpath==ik);
end
% infer transition matrix:
P = zeros(K);
t0=0;
for iseg=1:length(T)
    vpath_segment = vpath((t0+1):(t0+T(iseg)));
    for i1=1:K
        for i2=1:K
            P(i1,i2) = sum(vpath_segment(1:end-1)==i1 & vpath_segment(2:end)==i2);
        end
    end
    t0 = t0+T(iseg);
end
P = P./repmat(sum(P,2),1,K);
P_sim = cumsum(P,2);
% and simulate a new chain:
t0=0;
for iseg=1:length(T)
    z = rand(T(iseg),1);
    vpath_new(t0+1) = find(z(1)<cumsum(FO,2),1);
    for iT = 2:T(iseg)
        try % wierd syntax here just to catch very rare precision errors
            vpath_new(iT) = find(z(iT)<P_sim(vpath_new(iT-1),:),1);
        catch
            vpath_new(iT) = find(z(iT)<P_sim(vpath_new(iT-1),:),1);
        end
    end
    t0 = t0+T(iseg);
end
end