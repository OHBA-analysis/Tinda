function angleplot = circle_angles(bestseq)

% save each subject's rotational strength by this metric:
K=length(bestseq);
disttoplot_manual = zeros(K,1);
for i3=1:12
  disttoplot_manual(bestseq(i3)) = exp(sqrt(-1)*i3/K*2*pi);
end
angleplot = exp(sqrt(-1)*(angle(disttoplot_manual.')-angle(disttoplot_manual)));

