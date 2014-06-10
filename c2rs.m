function [x1, x2, x3] = c2rs(z)
% C2RS Cartesian complex coordinate to Riemann sphere projection.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.

theta = angle(z);
absz = abs(z);
phi = atan2(absz.^2 - 1, 2*abs(z));
phi(isinf(z)) = pi/2;
[x1, x2, x3] = sph2cart(theta, phi, ones(size(theta)));
