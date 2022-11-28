function rotational_momentum = compute_rotational_momentum(angleplot, FO_assym, whichstate)
tmp = angleplot.*FO_assym;
if exist('whichstate','var') && ~isempty(whichstate)
  % Note that we are counting each (i,j) double because for the rotational
  % momentum per state we take into account (i,j) and (j,i) for all j and one
  % particular i.
  tmp = squeeze(tmp(whichstate,:,:)) + squeeze(tmp(:,whichstate,:));
  rotational_momentum = imag(sum(tmp));
else
  rotational_momentum = imag(sum(sum(tmp)));
end
rotational_momentum = squeeze(rotational_momentum);