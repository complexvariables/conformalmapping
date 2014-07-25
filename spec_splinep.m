%% Specification by example -- periodic spline.
% This file best viewed via the command
%
%   open(publish('spec_splinep')).
%
clear

%%
% One way to create a spline object is to give a set of knots in the plane.
s = splinep([ ...
    0.5896 + 1.2486i; -0.1426 + 1.5954i; -0.9133 + 1.1561i
    -0.8465 + 0.3536i; -1.1116 - 0.2398i; -1.2695 - 0.9643i
    -0.5660 - 1.1075i; 0.2013 - 0.7552i; 0.8362 - 0.9634i
    1.5838 - 0.7013i; 1.3141 + 0.4008i; 0.8474 + 0.7291i ...
]);
plot(s)
hold on
plot(s.zpts, 'rd') % plot the knots

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
