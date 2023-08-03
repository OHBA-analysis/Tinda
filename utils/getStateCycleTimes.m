function cycletimes = getStateCycleTimes(vpath,T)

if size(vpath,2)>1
    % convert from Gamma to vpath:
    [~,vpath] = max(vpath,[],2);
end
num_fullcycles = 0;
time_tot = 0;
cycletimes = [];
for i=1:length(T)
    seq = vpath((sum(T(1:i-1))+1):sum(T(1:i)));
    firstchange = find(diff(seq)==1 | diff(seq)==(-max(seq)+1),1);
    if ~isempty(firstchange)
        fullcycles = find(seq(1:end-1)==seq(firstchange) & seq(2:end)==seq(firstchange+1));
        %num_fullcycles = num_fullcycles + length(fullcycles);
        %time_tot = time_tot + fullcycles(end)-fullcycles(1);
        cycletimes = [cycletimes;diff(fullcycles)];
    end
end

end