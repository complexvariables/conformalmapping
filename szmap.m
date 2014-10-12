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
    theKernel
    coefficients
    sopts
end

methods
    function f = szmap(range, varargin)
        if nargin
            if isa(range, 'closedcurve')
                range = region(range);
            end
            if ~issimplyconnected(range)
                error('CMT:InvalidArgument', ...
                    'Region must be simply connected.')
            end
            args = {unitdisk, range};
        else
            args = {};
        end

        f = f@conformalmap(args{:});
        opts = get(f, szset);
        
        if ~nargin
            return
        end

        if nargin > 2
            opts = set(opts, varargin{:});
        end

        boundary = outer(range);
        
        szargs = varargs(opts);
        S = szego(boundary, szargs{:});
        
        nF = opts.numFourierPts;
        t = invtheta(S, 2*pi*(0:nF-1)'/nF);
        c = flipud(fft(boundary(t))/nF);

        f.theKernel = S;
        f.coefficients = c;
        f.sopts = opts;
    end

    function g = ctranspose(f)
        d = domain(f);
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
        
        prefs = varargs(f.sopts);
        a = f.sopts.confCenter;
        g = szmap(cinvcurve(br, a), prefs{:}, 'confCenter', 0);
        func = @(z) 1./conj(applyMapPrivate(g, conj(1./z))) + a;

        g.theDomain = d;
        g.theRange = r;
        g.functionList = {func};
    end

    function S = kernel(f)
        S = f.theKernel;
    end
end

methods(Access=protected)
    function w = applyMap(f, z)
        if isanonymous(f)
            w = applyMap@conformalmap(f, z);
        else
            w = applyMapPrivate(f, z);
        end
    end
end

methods(Access=private)
    function w = applyMapPrivate(f, z)
        w = polyval(f.coefficients, z);
    end
end

end
