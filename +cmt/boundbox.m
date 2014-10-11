function box = boundbox(points)
% BOUNDBOX calculates the bounding box around points.
%
% box = boundbox(points) calculates a bounding box around a set of points
% in the complex plane, and returns coordinates in AXIS format.
%
% See also axis, plotbox.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.

box([1 3]) = min([real(points(:)) imag(points(:))], [], 1);
box([2 4]) = max([real(points(:)) imag(points(:))], [], 1);

end
