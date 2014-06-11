%% Specification by example -- periodic spline.
% This file best viewed via the command
%
%   open(publish('spec_mobius_grids')).
%
clear

%%
% One way to create a spline object is to give a set of knots in the plane.
s = splinep([ ...
    0.5896 + 1.2486i; -0.1426 + 1.5954i; -0.9133 + 1.1561i
    -0.5819 + 0.3160i; -0.8748 - 0.1156i; -1.5222 - 0.6936i
    -1.4220 - 1.4258i; -0.4971 - 0.8863i; 0.2274 - 0.4547i
    0.8362 - 0.9634i; 1.5838 - 0.7013i; 1.3141 + 0.4008i
    0.5742 + 0.4855i ...
]);
plot(s)
hold on
plot(zpts(s), 'rd') % plot the knots

%%
% If you are not into tediously entering knot points by hand, you may use
% the command (note the empty array as the only argument)
%
%   s = splinep([]);
%
% which will start an interactive figure in which you choose the knot
% points with the mouse button. Press the < enter > key after selecting the
% last point, or use the right mouse button to select the last point. You
% can then use the command
%
%   replicate(s)
%
% to print something you can copy and paste into a script. As was done with
% the given set of knots above. You should try it!


%%
% Splines, or rather their knots, may be subjected to arbitrary affine
% transformations.
clf
subplot(1,2,1)
plot(50*exp(1i*pi/2)*s)
subplot(1,2,2)
plot(s + 70 - 25i);


%%
% Splines may even be subjected to Mobius transformations. (You should
% probably be VERY careful with this.)
m = mobius([2, 2i, -2], [-1, -2i, 0]);
clf
plot(m(s))
