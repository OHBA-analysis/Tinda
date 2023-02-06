function idx = replay_scores2index(vpath, replayScores, upsample_factor, t_window, percentile)
if ~exist('vpath', 'var')
  vpath = [];
end
if ~exist('upsample_factor', 'var') || isempty(upsample_factor)
  upsample_factor = 2.5;
end
if ~exist('t_window', 'var') || isempty(t_window)
  t_window = 125;
end
if ~exist('percentile', 'var') || isempty(percentile)
  percentile = 1;
end

nSj = size(replayScores,1);

for iSj = 1:nSj
  if isempty(vpath)
    ntimes = 75000;
  else
    ntimes = length(vpath{iSj});
  end
  replaytimes = replayScores(iSj,:) > prctile(replayScores(iSj,:),100-percentile);
  
  replaytimes= [0,diff(replaytimes)]==1;
  %eliminate border points:
  replaytimes(1:t_window)=0;replaytimes([end-t_window]:end)=0;
  t_i = find(replaytimes);
  t_g = round(t_i*upsample_factor);
  
  if any(t_g>ntimes-t_window)
    t_g(t_g>ntimes-t_window) = [];
    t_i(t_i>ntimes-t_window) = [];
  end
  
  idx{iSj} = t_g;
end