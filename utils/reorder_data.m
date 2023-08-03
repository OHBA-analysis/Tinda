function data_new_order = reorder_data(orig_order, new_order, data_orig_order)
% this function takes in two vectors of with some ordering and a data
% matrix of any size. The data will be permuted along the relevant axis to
% move from the orig_order to the new_order.

sz = size(data_orig_order);
dim = find(sz==length(orig_order));

% permute to get the dimension that needs reordered first,
if dim~=1
  other_dims = setdiff(1:length(sz),dim);
  data_orig_order = permute(data_orig_order, [dim, other_dims]);
end


for i=1:length(new_order)
  mapping(i) = find(orig_order == new_order(i));
end

sz2 = size(data_orig_order);
data_new_order = reshape(data_orig_order(mapping,:), sz2);

% permute back to original dimensions
if dim~=1
  if dim==length(sz)
    data_new_order = permute(data_new_order, [2:dim 1]);
  else
    data_new_order = permute(data_new_order, [2:dim 1 dim+1:length(sz)]);
  end
end



