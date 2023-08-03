function [matched_id, F] = match_states(coh_template, coh_to_match)
% based on https://github.com/OHBA-analysis/osl-dynamics/blob/main/osl_dynamics/inference/modes.py#L205
% matching states based on correlations of coherence maps

n_matrices = size(coh_template,1);
F = zeros(n_matrices);
for i=1:n_matrices
    for j=1:n_matrices
        F(i,j) = corr(squash(triu(squeeze(coh_template(i,:,:)))), squash(triu(squeeze(coh_to_match(j,:,:)))));
    end
end

matched_id = matchpairs(F', 0, 'max');
matched_id = matched_id(:,1);