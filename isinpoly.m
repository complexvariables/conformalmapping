function in = isinpoly(z, vertex)
% ISINPOLY checks if point z is in polygon defined by vertex array.
%
%   Complex variable wrapper for builtin INPOLYGON.
%
%   See also INPOLYGON.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.

in = inpolygon(real(z),imag(z),real(vertex),imag(vertex));
