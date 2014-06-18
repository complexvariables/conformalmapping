classdef szmap < conformalmap
% SZMAP represents a Riemann map via the Szego kernel.
%
% f = szmap(range, a)
% f = szmap(range, a, opts)
% Creates a Riemann map from the interior of the unit disk to the interior of
% the region via the Szego kernel and the FFT with normalization f(0) = a and
% f'(0) > 0. Various options may be set via the opts parameter which must be a
% szset object.
%
% See also szego, szset.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

properties
  kernel_
  coefs_
  opts_
end

methods
  function f = szmap(range, a, opts)
    if ~nargin
      return
    end
    
    if isa(range, 'closedcurve')
      range = region(range);
    end
    if ~issimplyconnected(range)
      error('CMT:InvalidArgument', ...
            'Region must be simply connected.')
    end
    
    f = f@conformalmap(unitdisk, range);
    
    if nargin < 3 || isempty(opts)
      opts = szset;
    elseif ~isa(opts, 'szset')
      error('CMT:InvalidArgument', ...
            'Second argument must be a szset object.')
    end

    boundary = outer(range);
    S = szego(boundary, a, opts);
    nF = opts.nF;
    t = invtheta(S, 2*pi*(0:nF-1)'/nF);
    c = flipud(fft(boundary(t))/nF);
    
    f.kernel_ = S;
    f.coefs_ = c;
    f.opts_ = opts;
  end
  
  function g = ctranspose(f)
    d = domain(f);
    if ~issimplyconnected(d)
      error('CMT:NotDefined', ...
            'This operation only defined for simply connected domains.')
    end
    
    if isinterior(d)
      d = diskex(outer(d));
    else
      d = disk(inner(d));
    end
    
    r = range(f);
    if isinterior(r)
      br = outer(r);
      r = region(br, 'exteriorto');
    else
      br = inner(r);
      r = region(br);
    end
    
    stmp = szmap(br', 0, f.opts_);
    ftmp = @(z) 1./conj(apply_map(stmp, conj(1./z)));
    g = conformalmap(d, r, ftmp);
  end
end

methods(Access=protected)
  function w = apply_map(f, z)
    w = polyval(f.coefs_, z);
  end
end

end
