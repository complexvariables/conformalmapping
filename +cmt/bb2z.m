function z = bb2z(box)
%BB2Z axis bounding box to vertices.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.

z = complex(box([1, 2, 2, 1]), box([3, 3, 4, 4]));
z = z(:);

end
