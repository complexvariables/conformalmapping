function b = isdc(r)

b = ~isempty(r.innerboundary) & ~isempty(r.outerboundary);
 