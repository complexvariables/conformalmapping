function fill(d)

newplot
washold = ishold;

% Draw inner boundary
if ~isempty(d.innerboundary)
  hin = fill(d.innerboundary);
  zin = get(hin,'xdata') + 1i*get(hin,'ydata');
  hold on
  set(hin,'facecolor','none')
else
  zin = [];
end

if ~isempty(d.outerboundary)
  hout = fill(d.outerboundary);
  zout = get(hout,'xdata') + 1i*get(hout,'ydata');
  hold on
else
  axlim = axis;
  midpt = [mean(axlim(1:2)) mean(axlim(3:4))];
  width = [diff(axlim(1:2)) diff(axlim(3:4))];
  xbox = midpt(1) + width(1)*0.75*[-1 1];
  ybox = midpt(2) + width(2)*0.75*[-1 1];
  zout = xbox([1 1 2 2 1]) + 1i*ybox([2 1 1 2 2]);
end

if ~isempty(zin)
  set(hout,'facecolor','none')
  zall = [zout zin(end:-1:1)];
  settings = { [0.75 0.75 0.85],'edgecolor','b','linewidth',1.5, varargin{:} };
  h = fill(real(zall),imag(zall),settings{:});
end