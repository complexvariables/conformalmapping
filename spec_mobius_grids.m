%% Specification by example -- Mobius transformations.
% This file best viewed via the command
%
%   open(publish('spec_mobius_grids')).
%
clear


%%
% Create a basic Mobius transformation.
z3 = [1, 1i, -1];
w3 = [-1, -2i, 0];
m1 = mobius(z3, w3);

%%
% The Mobius transform is a conformal map.
if isa(m1, 'conformalmap')
  disp('m1 is a conformal map.')
else
  disp('m1 is not a conformal map.')
end

%%
% Plotting a conformal map gives an image in the range of a grid in the
% domain.
plot(m1)

%%
% Let's see that again in slow motion.
% (NOTE:
% Technically a Mobius map is an entire function, but for visualization
% convenience it uses the 3-vectors it was given on construction to define
% circles for its domain and range.)
gd = grid(domain(m1));
clf, hold on
plot(m1(gd), 'color', cmtplot.grey)
plot(range(m1), 'color', 'k')
axis(plotbox(range(m1)))
aspectequal
axis off


%%
% Another transformation for a composition example.
s3 = [0, 3, 1i];
m2 = mobius(w3, s3);

%%
% The |*| operator is the composition operator for conformal maps, with
% each map applied in the standard order. The following creates a new map
% which applies |m1|, then |m2|. Internally Mobius multiplies the matrices for
% |m1| and |m2| and creates a new Mobius object with the resulting matrix.
% This explains the discrepancy with the result of the example below.
m = m2*m1;
fprintf('composition comparison: %.g\n', norm(m2(m1(z3)).' - m(z3).'))

%%
% The base class has composition constructor behavior. Given arguments of the
% form (map1, map2, ..., mapn), it creates a new map via the composition
%    |mapn * ... * map2 * map1|.
% Internally |conformalmap| stores an array of the given maps and applies
% each in order.
mf = conformalmap(m1, m2);
fprintf('other comparison: %.g\n', norm(m2(m1(z3)).' - mf(z3).'))
