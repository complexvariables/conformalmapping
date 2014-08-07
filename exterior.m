function R = exterior(C)
% EXTERIOR creates an exterior region bounded by curve.
%
% R = exterior(C)
% Creates exterior region R bounded by closed curve C.
%
% See also closedcurve, region.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.

if ~isa(C, 'closedcurve')
    error('CMT:InvalidArgument', ...
        'Function argument must be a closed curve.')
end

R = region(C, 'exteriorto');
