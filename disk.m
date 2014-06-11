classdef disk < region
% DISK is a region bounded by a circle.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% (Re)written by Everett Kropf, 2014,
% adapted from an idea by Toby Driscoll, 20??.

methods
  function D = disk(center, radius)
    badargs = false;
    switch nargin
      case 0
        C = [];
        
      case 1
        if isa(center, 'disk')
          C = center.outerboundary_;
        elseif isa(center, 'double') && numel(center) == 3
          C = circle(center);
        elseif isa(center, 'circle') && ~isinf(center)
          C = center;
        else
          badargs = true;
        end
        
      case 2
        if isa(center, 'double') && isa(radius, 'double') ...
            && numel(center) == 1 && numel(radius) == 1
          z3 = center + radius*exp(2i*pi*[0, 0.25, 0.5]);
          C = circle(z3);
        else
          badargs = true;
        end
        
      otherwise
        badargs = true;
    end
    if badargs
      error('CMT:InvalidArgument', ...
            'Expected 3 points or a center and radius.')
    end
    
    if isempty(C)
      supargs = {};
    else
      supargs = {C};
    end
    D = D@region(supargs{:});
  end
  
  function gd = grid(D, radial, circular)
    if nargin < 2 || isempty(radial)
      nrad = 20;
    else
      nrad = radial;
    end
    if nargin < 3 || isempty(circular)
      ncirc = 5;
    else
      ncirc = circular;
    end
    
    npt = 200;
    c = center(outer(D));
    r = radius(outer(D));
        
    curves = cell(nrad + ncirc, 1);
    zg = (1:npt)'/(npt+1);
    for k = 1:nrad
      curves{k} = c + r*exp(2i*pi*(k-1)/nrad)*zg;
    end
    zg = exp(2i*pi*(0:npt-1)'/(npt-1));
    for k = 1:ncirc
      curves{nrad + k} = c + r*k/(ncirc+1)*zg;
    end
    
    gd = gridcurves(curves);
  end
end

end
