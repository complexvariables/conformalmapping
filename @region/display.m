function display(d)

if isempty(d)
  fprintf('\nempty region\n\n')
  return
end

fprintf('\nregion')

if ~isempty(d.outerboundary)
  fprintf(' interior to:\n')
  display(d.outerboundary)
  if ~isempty(d.innerboundary)
    fprintf('\n and')
  end
end

if ~isempty(d.innerboundary)
  fprintf(' exterior to:\n')
  display(d.innerboundary)
end

