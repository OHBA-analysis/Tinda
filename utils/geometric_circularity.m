function [circularity, pval, circularity_perm, permuted_graphs, fig] = geometric_circularity(mean_direction, sigpoints, nperm, doplot)
% this function provides a measure of how geometrical circular a directed
% graph is. It is based on the directed (i.e., clockwise vs. counter
% clockwise) distance between connected nodes. For perfect circles, the
% directed distance (clockwise or counterclockwise) is small, because the
% connections lie close to each other on the unit circle, in the direction
% that is the same as the mean circle direction. For example, if the mean
% direction is counterclockwise, going from pi/4 to pi/2 is a small
% distance (pi/2), whereas going from pi/2 to pi/4 is a large distance (7/4
% pi).
% In order to test whether the directed graph is more circular then
% expected by chance (given the number of directed edges), we compare it
% with the circularity measure after permuting edge (directions). We still
% have to figure out what the most fair permutation scheme is. Naively, we
% expect a random graph to have average jumps of pi.
%
% INPUT
% mean_direction:   N x N matrix of (directed) edge weight
% sigpoints:        N x N matrix containing ones for statistically
%                   significant edges
% nperm:            double/int, number of permutations
% doplot:           boolean, whether to plot empirical and permuted 
%                   circularity
%
% OUTPUT
% circularity:        double - measure between 0 (low circularity) and 1 
%                     (high circularity), defined by a normalized measure
%                     of empirical circularity over mean over permutations
% pval:               double, p-value for empirical circularity vs. 
%                     permutations
% circularity_perm:   N x 1 vector containing permuted circularity scores
% permuted_graphs:    cell(N,1) containing the permuted graphs (K x K)
% fig:                figure object
%
%
% Copyright Mats van Es, University of Oxford, August 2022.

if nargin<3 || isempty(nperm)
  nperm = 1000;
end
if nargin<4 || isempty(doplot)
  doplot = 1;
end

% create some (geometrical) variables
K = size(mean_direction,1);
phase = zeros(K,1);
for i=1:K
  phase(i) = exp(sqrt(-1)*(i+2)/K*2*pi);
end
offdiag=find(~eye(K));
nsigpoints = sum(sigpoints(:));

sigpoints_sign = sign(mean_direction).*sigpoints; % get the signed
% significant connections
clockdir = sign(mean(sign(mean_direction(:))));% determine whether the mean
% circle direction is clockwise (1) or counter clockwise (-1)

% define directed distance
dist = circ_dist2(angle(phase), angle(phase)); % distances between points
% on unit circle (bidirectional) NOTE: relies on circular statistics
% toolbox

dist(dist<0) = 2*pi+dist(dist<0); % transform counterclockwise distance
% into clockwise distance (i.e. small counterclockwise distance means large
% clockwise distance) in case the mean direction is counter clockwise, we
% want the reverse to be true:
if clockdir==-1
  dist = dist';
end

dist = dist - diag(diag(dist)); % set 2*pi on diagonal to zero

% create random permutations for offdiagonal elements
for k=1:nperm
  P(k,:) = randperm(K.^2-K); % randomly permuting arrows on the circle plot
  P2(k,:) = [1,randperm(11)+1]; % this is Cam's approach
  P3(k,:) = rand(1, (K.^2-K)/2); % randomly swap signs between i,j and j,i
end
P3=round(P3);


mean_phasejump_perm = zeros(nperm,1);
for k=1:nperm+1
  sigpoints_sign_tmp = sigpoints_sign; % reformat the signed sigpoints
  % according to the overall pattern
  if k>1
    % permute offdiagonal elements
    if 1
      sigpoints_sign_tmp(offdiag) = sigpoints_sign_tmp(offdiag(P(k-1,:)));
      
    elseif 1
      % this is Cam's approach
      sigpoints_sign_tmp = sigpoints_sign(P2(k-1,:), P2(k-1,:));
      
    elseif 1
      % another way to do it is by exchanging the sign between i,j and j,i
      % perhaps?
      x0 = tril(ones(K),-1);
      cnt=1;
      for ik1=1:K
        for ik2=setdiff(1:K,ik1)
          if x0(ik1,ik2)
            if P3(k-1, cnt)
              sigpoints_sign_tmp(ik1,ik2) = sigpoints_sign(ik2, ik1);
              sigpoints_sign_tmp(ik2,ik1) = sigpoints_sign(ik1, ik2);
            end
            cnt=cnt+1;
          end
        end
      end
    end
  end
  % seperately get the distances for the sigpoints that are positive (i->j)
  % and negative (j->i). A large distance for i->j means small distance for
  % j->i
  permuted_graphs{k,1} = sigpoints_sign_tmp;
  tmp1 = sigpoints_sign_tmp; tmp2 = sigpoints_sign_tmp;
  tmp1(tmp1<0)=0; tmp2(tmp2>0)=0; tmp2(tmp2<0)=1;
  
  % sum the distances
  d = (sum(sum(tmp1.*dist)) + sum(sum(tmp2.*transpose(dist))))./nsigpoints;
  if k==1
    mean_phasejump = d;
  else
    mean_phasejump_perm(k-1,1) = d;
  end
end

circularity = compute_circularity(mean_phasejump, mean_phasejump_perm, K);
pval = max([(nperm - sum(mean_phasejump < mean_phasejump_perm)) / nperm, 1 / nperm]);
circularity_perm = compute_circularity(mean_phasejump_perm,mean_phasejump_perm,K);

if doplot
  fig=figure; 
  
  % plot the cycle with the smallest phasejump in the permutation
  subplot(2,3,1)
  cyclicalstateplot(1:K,mean_direction, sigpoints, [], false);
  title(sprintf('empirical cycle plot \n circularity = %s \n', num2str(round(circularity, 2))))
  [~,ix] = min(mean_phasejump_perm);
  
  
  subplot(2,3,4)
  cyclicalstateplot(1:K,permuted_graphs{ix}, abs(permuted_graphs{ix}), [], false);
  title(sprintf('max permutation  circularity \n circularity = %s \n', num2str(round(max(circularity_perm),2))))
  
  subplot(2,3,[2,5])
  histogram(mean_phasejump_perm)
  hold on, vline(mean_phasejump), title(sprintf('empirical mean phase \n jump vs. permutations'))
  xlabel('Radians'), ylabel('count'), xlim([0 2*pi])
  
  subplot(2,3,[3,6])
  histogram(circularity_perm)
  hold on, vline(circularity), title(sprintf('empirical circularity \n vs. permutations'))
  xlabel('Radians'), ylabel('count'), xlim([0, 1.1])
else
  fig = false;
end

function circularity = compute_circularity(phasejump,phasejump_perm,K)
  circularity = (abs(circ_dist(pi, phasejump))+2*pi/K)/(2*pi-abs(circ_mean(phasejump_perm)));
end

end

