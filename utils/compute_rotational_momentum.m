function rotational_momentum = compute_rotational_momentum(angleplot, FO_assym, whichstate)
tmp = angleplot.*FO_assym;
nSj = size(FO_assym,3);
if exist('whichstate','var') && ~isempty(whichstate)
  % Note that we are counting each (i,j) double because for the rotational
  % momentum per state we take into account (i,j) and (j,i) for all j and one
  % particular i.
  tmp = squeeze(tmp(whichstate,:,:)) + squeeze(tmp(:,whichstate,:));
  rotational_momentum = imag(nansum(tmp));
else
  rotational_momentum = imag(nansum(nansum(tmp)));
end
rotational_momentum = squeeze(rotational_momentum);
if size(rotational_momentum,2) > size(rotational_momentum,1) && size(rotational_momentum,1)==1
  rotational_momentum = rotational_momentum';
end
rotational_momentum = -rotational_momentum; % positive rotational momentum should indicate clockwise cycle