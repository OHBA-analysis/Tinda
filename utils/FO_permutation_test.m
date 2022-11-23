function [pvals, stat] = FO_permutation_test(FO, K, nSj)

dat1=[];
dat1.dimord = 'rpt_chan_time';
dat1.label{1} = 'FO asym';

dat1.time=1:(K^2-K);
dat2=dat1;

cfg=[];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.design = [ones(1,nSj), 2*ones(1,nSj); 1:nSj, 1:nSj];
cfg.ivar = 1;
cfg.uvar = 2;
cfg.numrandomization = 100000;
cfg.tail = 0;
cfg.correcttail = 'prob';
pvals = zeros(K,K);

tmp = permute(squeeze(FO(:,:,1,:)), [3,1,2]);
dat1.trial(:,1,:) = tmp(:, ~eye(K));

tmp = permute(squeeze(FO(:,:,2,:)), [3,1,2]);
dat2.trial(:,1,:) = tmp(:, ~eye(K));

stat = ft_timelockstatistics(cfg, dat1, dat2);
tmp = ones(K);
tmp(~eye(K)) = stat.prob;
pvals(:,:) = tmp;
