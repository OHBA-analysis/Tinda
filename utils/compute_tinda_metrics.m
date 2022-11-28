function metrics = compute_tinda_metrics(config, bestseq, angleplot, FO_intervals, sigpoints, color_scheme, doplot)
if ~exist('doplot')
  doplot = true;
end
K = size(FO_intervals,1);
metrics = struct;

metrics.mean_direction = squeeze(mean(FO_intervals(:,:,1,:)-FO_intervals(:,:,2,:),4));
metrics.mean_assym = squeeze(mean((FO_intervals(:,:,1,:)-FO_intervals(:,:,2,:))./mean(FO_intervals,3),4));

metrics.FO_assym = squeeze((FO_intervals(:,:,1,:)-FO_intervals(:,:,2,:))./mean(FO_intervals,3));


% Rotational momentum is the rotational strength of the unthresholded FO
% asymmetry
rotational_momentum = compute_rotational_momentum(angleplot, metrics.FO_assym);
metrics.rotational_momentum = rotational_momentum;
metrics.max_theoretical_rotational_momentum = compute_rotational_momentum(angleplot, sign(imag(angleplot)));

% also get the metric for each state
for k=1:K
  metrics.rotational_momentum_perstate(:,k) = compute_rotational_momentum(angleplot, metrics.FO_assym, k);
end


% TIDA is the mean(abs(FO_assym))
metrics.TIDA = nanmean(reshape(abs(metrics.FO_assym), K*K,[]))';
% also compute this per state
for k=1:K
  metrics.TIDA_perstate(:,k) = nanmean(abs([squeeze(metrics.FO_assym(:,k,:));squeeze(metrics.FO_assym(k,:,:))]));
end


% Circularity is the normalised, thresholded clockwise distance (e.g. done
% on the circle plot
[circularity, circle_pval, ~, ~, fig] = geometric_circularity(bestseq, metrics.mean_direction, sigpoints,[],[],doplot,color_scheme);
metrics.circularity = circularity;
metrics.circularity_pval = circle_pval;
if doplot
  set_font(10, {'title', 'label'})
  save_figure([config.figdir,'2supp_CyclicalpatternVsPermutations']);
end

% get the measure per subject
tmp = metrics.FO_assym; 
tmp(isnan(tmp))=0;
for k=1:config.nSj
  [metrics.circularity_subject(k,1), metrics.circularity_pval_subject(k,1), ~, ~, ~] = geometric_circularity(bestseq, tmp(:, :,k), sigpoints,[],[],0, color_scheme);
end