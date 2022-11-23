function rotational_momentum = compute_rotational_momentum(angleplot, FO_assym, whichstate)
tmp = angleplot.*FO_assym;
if exist('whichstate','var') && ~isempty(whichstate) 
  tmp = squeeze(tmp(whichstate,:,:)) + squeeze(tmp(:,whichstate,:));
  rotational_momentum = imag(sum(tmp));
else
  rotational_momentum = imag(sum(sum(tmp)));
end
rotational_momentum = squeeze(rotational_momentum);