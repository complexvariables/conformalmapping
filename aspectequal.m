function aspectequal(ah)
% ASPECTEQUAL sets data aspect ratio to equal on all dimensions.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

if ~nargin || isempty(ah)
  ah = gca;
end

set(ah, 'dataaspectratio', [1, 1, 1])
