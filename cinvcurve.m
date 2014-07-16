classdef cinvcurve < closedcurve
% CINVCURVE class represents the conjugate inverse of a curve.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

properties
  curve_
end

methods
  function C = cinvcurve(curve)
    if ~nargin
      return
    end
    
    if ~isa(curve, 'closedcurve')
      error('CMT:InvalidArgument', 'Expected a closedcurve object.')
    end
    
    C.curve_ = curve;
  end
  
  function disp(C)
    fprintf('conjugate inverse of curve:\n')
    disp(C.curve_)
  end
  
  function z = point(C, t)
    z = 1./conj(point(C.curve_, t));
  end
  
  function zt = tangent(C, t)
    zt = -conj(tangent(C.curve_, t)./point(C.curve_, t).^2);
  end
end

end