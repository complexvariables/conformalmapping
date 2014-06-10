classdef corner
% CORNER class represents corners for closed curves.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

properties
  param
  point
  alpha
end

properties(Access=private)
  % Should probably use metadata.
  prop_list = {'param', 'point', 'alpha'}
end

methods
  function C = corner(cinfo)
    if ~nargin
      return
    end
    
    if isa(cinfo, 'corner')
      C = cinfo;
      return
    end
    
    % Allow cell and structure forms for construction.
    if isa(cinfo, 'struct')
      if ~all(isfield(cinfo, C.prop_list))
        error('Corner information was supplied incorrectly.')
      end
      for prop = C.prop_list
        C.(prop) = cinfo.(prop);
      end
    elseif isa(cinfo, 'cell')
      % Assume order of cell entries.
      for k = 1:numel(C.prop_list)
        C.(C.prop_list{k}) = num2cell(cinfo{k});
      end
    else
    end
  end
end

end