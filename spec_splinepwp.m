%% Specification by example -- periodic piecewise spline.
% This file best viewed via the command
%
%   open(publish('spec_splinepwp')).
%
clear


%%
% Some sample knots to get us started. We'll show a regular periodic
% spline with its knots numbered, and then add some corners.

knots = [1, 1i, -1, -1i];
G = splinep(knots);

plot(G), hold on
plot(knots, 'b.', 'markersize', 18)
for k = 1:numel(knots)
  text(real(knots(k)), imag(knots(k)), ['   ' int2str(k)])
end
axis off, hold off


%%
% Corners are easy to specify.

corners = [1, 2];
G = splinepwp(knots, corners);

plot(G), hold on
plot(knots, 'b.', 'markersize', 18)
for k = 1:numel(knots)
  text(real(knots(k)), imag(knots(k)), ['   ' int2str(k)])
end
axis off, hold off


%%
% Note that |splinepwp| requires a corner be the first knot, so the
% constructor will rearrange the knots as needed to make this happen.

corners = [3, 4];
G = splinepwp(knots, corners);

plot(G), hold on
plot(knots, 'b.', 'markersize', 18)
for k = 1:numel(G.knots)
  text(real(G.knots(k)), imag(G.knots(k)), ['   ' int2str(k)])
end
axis off, hold off


%%
% Interestingly, if the region isn't too "strange", the Szego kernel will
% give a map with boundary corners.

f = szmap(G, 0);
plot(f)
