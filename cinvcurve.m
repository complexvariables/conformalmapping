classdef cinvcurve < closedcurve
% CINVCURVE class represents the conjugate inverse of a curve.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

properties
  theCurve
  confCenter = 0
end

methods
  function C = cinvcurve(curve, confCenter)
    if ~nargin
      return
    end
    
    if ~isa(curve, 'closedcurve')
      error('CMT:InvalidArgument', 'Expected a closedcurve object.')
    end
    if nargin > 1
        if numel(confCenter) == 1 && isfinite(confCenter)
            C.confCenter = confCenter;
        else
            error('CMT:InvalidArgument', ...
                'Expected a finite scalar for a conformal center.')
        end
    end
    
    C.theCurve = curve;
  end
  
  function disp(C)
    fprintf('conjugate inverse of curve:\n')
    disp(C.theCurve)
  end
  
  function z = point(C, t)
    z = 1./conj(point(C.theCurve, t) - C.confCenter);
  end
  
  function zt = tangent(C, t)
    zt = -conj(tangent(C.theCurve, t))./...
        conj(point(C.theCurve, t) - C.confCenter).^2;
  end
end

end