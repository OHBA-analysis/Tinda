function set_sqrt_ax(f, f1, f2, x_or_y)
if ~exist('f1', 'var'), f1 = 1; end
if ~exist('f2', 'var'), f2 = floor(sqrt(f(end))); end
if ~exist('x_or_y', 'var'), x_or_y = 'x'; end

xt = (f1 : f2).^2;
sqrtf = sqrt(f);

for ik=1:length(xt)
  xtl{ik} = num2str(xt(ik));
  xt(ik) = nearest(f, xt(ik));
end
xt = sqrtf(xt);
xticks(xt)
if strcmp(x_or_y, 'x')
  set(gca, 'XTickLabel', xtl)
elseif strcmp(x_or_y, 'y')
  set(gca, 'YTickLabel', xtl)
end