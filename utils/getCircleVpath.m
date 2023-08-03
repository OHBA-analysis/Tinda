function vpathcircle = getCircleVpath(vpath, ordering)

disttoplot_manual = zeros(12,1);
for i=1:12
    temp = exp(sqrt(-1)*(i+2)/12*2*pi);
    disttoplot_manual(ordering(i),1) = temp;
end

vpathcircle = vpath;
for k=1:max(vpath)
    vpathcircle(vpath==k) = disttoplot_manual(k,:);
end