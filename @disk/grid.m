function g = grid(d,varargin)

[theta,r] = optargs({12,9},varargin);

% Convert positive integers into vectors of values.
if (length(theta)==1) & (theta>0) & (round(theta)==theta)
  theta = linspace(0,2*pi,theta+1);  theta(end) = [];
end
if (length(r)==1) & (r>0) & (round(r)==r)
  r = linspace(0,d.radius,r+2);
  r([1 end]) = [];
end

curve = {};

% Define the radii.
for j = 1:length(theta)
  curve{j} = @(t) d.center + d.radius*t*exp(1i*theta(j));
end

% Define the circles.
for j = 1:length(r)
  curve{end+1} = @(t) d.center + r(j)*exp(2i*pi*t);
end

g = gridcurves(d,curve{:});
