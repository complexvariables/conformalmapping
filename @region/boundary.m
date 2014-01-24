function b = boundary(r,select)

if nargin > 1 & isequal(select,'list')
  b = {r.outerboundary r.innerboundary};
elseif isempty(r.innerboundary)
  b = r.outerboundary;
elseif isempty(r.outerboundary)
  b = r.innerboundary;
else
  if nargin > 1
    if isequal(select,'inner')
      b = r.innerboundary;
    elseif isequal(select,'outer')
      b = r.outerboundary;
    else
      b.inner = r.innerboundary;
      b.outer = r.outerboundary;
    end
  end
end
