function metrics = compute_tinda_metrics(config, bestseq, FO_intervals, sigpoints, color_scheme, doplot, compute_circularity)
if ~exist('doplot')
  doplot = true;
end
if ~exist('compute_circularity')
  compute_circularity = false;
end
angleplot = circle_angles(bestseq);

K = size(FO_intervals,1);
metrics = struct;

metrics.mean_direction = squeeze(mean(FO_intervals(:,:,1,:)-FO_intervals(:,:,2,:),4));
metrics.mean_assym = squeeze(mean((FO_intervals(:,:,1,:)-FO_intervals(:,:,2,:))./mean(FO_intervals,3),4));
metrics.FO_assym = squeeze((FO_intervals(:,:,1,:)-FO_intervals(:,:,2,:))./mean(FO_intervals,3));
% also get measures for how well the fit is of subject with average pattern
tmp = metrics.FO_assym;
tmp(isnan(tmp))=0;
for k=1:size(tmp,3)
  metrics.FO_assym_subject_fit(k,1) = corr2(tmp(:,:,k), mean(tmp,3));
  metrics.FO_assym_rv_coef(k,1) = rv_coefficient(tmp(:,:,k), mean(tmp,3));
end

tmp = squeeze(mean((FO_intervals(:,:,1,:)./FO_intervals(:,:,2,:)),4));
metrics.P = tmp./nansum(tmp,2); % transition probability based on asymmetries



% Rotational momentum is the rotational strength of the unthresholded FO
% asymmetry
rotational_momentum = compute_rotational_momentum(angleplot, metrics.FO_assym);
metrics.rotational_momentum = rotational_momentum;
metrics.max_theoretical_rotational_momentum = abs(compute_rotational_momentum(angleplot, sign(imag(angleplot))));

% also get the metric for each state
% Note that we are counting each (i,j) double because for the rotational
% momentum per state we take into account (i,j) and (j,i) for all j and one
% particular i.
for k=1:K
  metrics.rotational_momentum_perstate(:,k) = compute_rotational_momentum(angleplot, metrics.FO_assym, k);
end
metrics.max_theoretical_rotational_momentum_perstate = abs(2*metrics.max_theoretical_rotational_momentum/K);



% TIDA is the mean(abs(FO_assym))
metrics.TIDA = nanmean(reshape(abs(metrics.FO_assym), K*K,[]))';
% also compute this per state
for k=1:K
  metrics.TIDA_perstate(:,k) = nanmean(abs([squeeze(metrics.FO_assym(:,k,:));squeeze(metrics.FO_assym(k,:,:))]));
end


% Circularity is the normalised, thresholded clockwise distance (e.g. done
% on the circle plot
if compute_circularity
  [circularity, circle_pval, ~, ~, fig] = geometric_circularity(bestseq, metrics.mean_direction, sigpoints,[],[],doplot,color_scheme);
  metrics.circularity = circularity;
  metrics.circularity_pval = circle_pval;
  if doplot
    set_font(10, {'title', 'label'})
    save_figure([config.figdir,'figure_supp_tinda_metrics/', '2supp_CyclicalpatternVsPermutations'],[],false);
  end
  
  % get the measure per subject
  tmp = metrics.FO_assym;
  tmp(isnan(tmp))=0;
  for k=1:config.nSj
    [metrics.circularity_subject(k,1), metrics.circularity_pval_subject(k,1), ~, ~, ~] = geometric_circularity(bestseq, tmp(:, :,k), sigpoints,[],[],0, color_scheme);
  end
end