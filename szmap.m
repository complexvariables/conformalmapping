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
    function f = szmap(range, a, varargin)
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
        S = szego(boundary, a, szargs{:});
        
        nF = opts.numFourierPts;
        t = invtheta(S, 2*pi*(0:nF-1)'/nF);
        c = flipud(fft(boundary(t))/nF);

        f.kernel_ = S;
        f.coefs_ = c;
        f.opts_ = opts;
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
        
        prefs = varargs(f.opts_);
        g = szmap(br', 0, prefs{:});
        func = @(z) 1./conj(apply_map_(g, conj(1./z)));

        g.domain_ = d;
        g.range_ = r;
        g.function_list_ = {func};
    end

    function S = kernel(f)
        S = f.kernel_;
    end
end

methods(Access=protected)
    function w = apply_map(f, z)
        if isanonymous(f)
            w = apply_map@conformalmap(f, z);
        else
            w = apply_map_(f, z);
        end
    end
end

methods(Access=private)
    function w = apply_map_(f, z)
        w = polyval(f.coefs_, z);
    end
end

end
