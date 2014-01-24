function plot(g,varargin)

if isempty(g.curve), return, end

newplot
washold = ishold;

% Two choices: functions that define the curves, or points.
% Find the bounding box from the data.
axlim = [Inf -Inf Inf -Inf];
for j = 1:length(g.curve)
  if isa(g.curve{j},'function_handle')
    [t,z{j}] = adaptplot(g.curve{j},[100*eps 1-100*eps]);
  else
    z{j} = g.curve{j};
  end
  axlim(1) = min( axlim(1), min(real(z{j})) );
  axlim(2) = max( axlim(2), max(real(z{j})) );
  axlim(3) = min( axlim(3), min(imag(z{j})) );
  axlim(4) = max( axlim(4), max(imag(z{j})) );
end

axis(axlim)
hold on

plot(g.region)
for j = 1:length(z)
  plot(real(z{j}),imag(z{j}),varargin{:});
end

if ~washold
  hold off
end
