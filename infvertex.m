function v = infvertex(zin, zout)
%INFVERTEX Create a representation of a vertex at infinity.
%   INFVERTEX(ZIN,ZOUT) creates an object that represents a vertex at
%   infinity, with an incoming straight side parallel to ZIN and an
%   outgoing straight side parallel to ZOUT.
%
%   For usage examples see POLYGON.
%
%   See also POLYGON/POLYGON, HOMOG.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.

v = homog([zin, zout], 0);
