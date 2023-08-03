
psd_wb_norm = psd_wb-mean(psd_wb,2);
coh_wb_norm = coh_wb-mean(coh_wb,2);
for i=1:size(psd_wb,1)
  for j=1:12
    for k=1:12
      r_coh_norm(i,j,k) = corr(squash(coh_wb_norm(i,j,:)), squash(coh_wb_norm(i,k,:)));
    end
  end
end
for i=1:size(psd_wb,1)
  for j=1:12
    for k=1:12
      r_coh(i,j,k) = corr(squash(coh_wb(i,j,:)), squash(coh_wb(i,k,:)));
    end
  end
end
for i=1:size(psd_wb,1)
  for j=1:12
    for k=1:12
      r_pow(i,j,k) = corr(squash(psd_wb(i,j,:)), squash(psd_wb(i,k,:)));
    end
  end
end
for i=1:size(psd_wb,1)
  for j=1:12
    for k=1:12
      r_pow_norm(i,j,k) = corr(squash(psd_wb_norm(i,j,:)), squash(psd_wb_norm(i,k,:)));
    end
  end
end



q=hmm_1stlevel.cycle_metrics.mean_assym
q(find(eye(12)))=0;

for k=1:size(psd_wb,1)
  R{1}(k,:) = corr(q(:), squeeze(r_pow(k,:))');
  R{2}(k,:) = corr(q(:), squeeze(r_coh(k,:))');
  R{3}(k,:) = corr(q(:), squeeze(r_pow_norm(k,:))');
  R{4}(k,:) = corr(q(:), squeeze(r_coh_norm(k,:))');
end

%%
r_coh_norm_mean = squeeze(mean(r_coh_norm));

% r_coh_norm_mean = q2.*abs(hmm_1stlevel.cycle_metrics.mean_assym);
for k=1:12
if k<12
Y1(k) = r_coh_norm_mean(k,k+1);
Y2(k) = r_coh_norm_mean(bestseq(k), bestseq(k+1));
else
Y1(k) = r_coh_norm_mean(12,1);
Y2(k) = r_coh_norm_mean(bestseq(12), bestseq(1));
end
end
sum_correlation_default = sum(Y1)
sum_correlation_bestseq = sum(Y2)


nperm=10000;
for iperm=1:nperm
bestseq_perm = randperm(12);
for k=1:12
if k<12
Y_perm(k,iperm) = r_coh_norm_mean(bestseq_perm(k), bestseq_perm(k+1));
else
Y_perm(k,iperm) = r_coh_norm_mean(bestseq_perm(12), bestseq_perm(1));
end
end
end
sum(sum(Y_perm)>sum_correlation_bestseq)/nperm

