function box = plotbox(points, scale)
% PLOTBOX returns padded axis coordinates around points.
%
% box = plotbox(points) calculates a 1-by-4 array to pass to AXIS which
% sets a padded square box around the given points.
%
% box = plotbox(points, scale) allows setting the padding scale, which
% defaults to 1.2 times the largest axis of the bounding box.
%
% See also axis, boundbox.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

if nargin < 2 || isempty(scale)
  scale = 1.2;
end

box = boundbox(points);

dbox = scale/2*max(diff(box(1:2)), diff(box(3:4)))*[-1 1];

box(1:2) = mean(box(1:2)) + dbox;
box(3:4) = mean(box(3:4)) + dbox;