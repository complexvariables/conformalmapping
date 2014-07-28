classdef splinepwp < closedcurve
% SPLINEPWP is a periodic piecewise spline with corners.
%
% G = splinepwp(knots, corners)
% The corners array is a set of indicies of the knots array which dictates
% which knots in the list will be corners. The array corners should be
% strictly monotonically increasing with
%   1 <= corners(1) and corners(end) <= numel(knots).
% If corners(1) ~= 1, then the knots and corners arrays will be rearranged
% so that this is the case.
%
% Currently the corner angles are determined strictly by the angles of
% neighboring knots, by pretending we are making a polygon with the knots
% as vertices.
%
% This class is mainly for testing corner cases on boundaries that aren't
% polygons. For instance it is completely arbitrary that the value of the
% tangent at a corner has the direction of the outgoing tangent, as it
% were. In fact to create the splines between corners, the endpoint tangent
% vectors arbitrarily have a magnitude of 5, which was selected because it
% gives the "nicest, natural" distribution of points for an evenly spaced
% discrete parameterization.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

properties(SetAccess=protected)
    knots
    corners
    alpha
    tknots
    cornerTangent
end

properties(Access=protected)
    ppx
    ppy
end

methods
    function G = splinepwp(knots, corners)
        % Constructor.

        if ~nargin
            return
        end

        % Check input.
        corners = corners(:);
        knots = knots(:);
        if knots(1) == knots(end)
            knots = knots(1:end-1);
        end
        if any(diff(corners) <= 0)
            error('CMT:BadArguemnt', ...
                'The corners array must be strictly monotonically increasing.')
        end
        try
            knots(corners);
        catch err
            if strcmp(err.identifier, 'MATLAB:badsubscript')
                error('CMT:BadArgument', ...
                    'The corners array must be a valid set of indices for knots.')
            else
                rethrow(err)
            end
        end

        % Constructor work.
        G.knots = knots;
        G.corners = corners;
        G = calculateSplines(G);
    end

    function tk = breaks(G)
        tk = [G.tknots(G.corners); 1];
    end

    function z = point(G, t)
        % Point on curve at t.

        z = paramEval(G, t, 0);
    end

    function zt = tangent(G, t)
        % Tangent to curve at t.

        zt = paramEval(G, t, 1);
    end
end

methods(Access=protected)
    function G = calculateSplines(G)
        v = G.knots;
        cv = G.corners;
        if cv(1) ~= 1
            % Shift corners and knots.
            shift = cv(1) - 1;
            v = circshift(v, -shift);
            cv = cv - shift;
        end

        % Pseudo arc length parameterization via chordal arclengths.
        nv = numel(v);
        v = [v; v(1)];
        dL = abs(diff(v));
        t = [0; cumsum(dL)];
        t = t/t(end);

        % Corner angles determined by neighboring knots.
        pv = mod(cv - 2, nv) + 1;
        fv = mod(cv, nv) + 1;
        alpha_ = mod(angle((v(pv) - v(cv))./(v(fv) - v(cv)))/pi, 2);

        % Outgoing/incoming tangent vectors.
        vtan = [ v(fv) - v(cv), v(cv) - v(pv) ];
        mu = 5; % (arbitrary) tangent vector magnitude
        vtan = mu*vtan./abs(vtan);

        % Compute piecewise polynomial splines for segments between corners.
        ncv = numel(cv);
        ppx_(ncv, 3) = mkpp([1 1], [1 1]);
        ppy_ = ppx_;
        for k = 1:ncv
            if k + 1 <= ncv
                vdx = cv(k):cv(k+1);
            else
                vdx = cv(k):nv+1;
            end
            vex = [vtan(k,1); v(vdx); vtan(mod(k, ncv)+1,2)];
            tk = t(vdx);

            % Piecewise polynomials and derivatives.
            ppx_(k,1) = spline(tk, real(vex));
            ppx_(k,2:3) = splinepwp.ppderivs(ppx_(k,1));
            ppy_(k,1) = spline(tk, imag(vex));
            ppy_(k,2:3) = splinepwp.ppderivs(ppy_(k,1));
        end

        G.knots = v(1:end-1);
        G.corners = cv;
        G.alpha = alpha_;
        G.ppx = ppx_;
        G.ppy = ppy_;
        G.tknots = t;
        G.cornerTangent = vtan;
    end

    function z = paramEval(G, t, d)
        % Evaluate curve given parameter using derivative level d.

        if nargin < 3
            % No derivative.
            d = 0;
        end
        d = d + 1;
        t = modparam(G, t);
        z = nan(size(t));

        brks = breaks(G);
        for k = 1:numel(G.corners)
            tL = brks(k) <= t & t < brks(k+1);
            z(tL) = complex(ppval(G.ppx(k,d), t(tL)), ppval(G.ppy(k,d), t(tL)));
        end
    end
end

methods(Static, Access=private)
    function ppd = ppderivs(pp)
        % First and second derivatives of piecewise polynomial pp.
        coef = pp.coefs;
        ppd = [ mkpp(pp.breaks, [3*coef(:,1), 2*coef(:,2), coef(:,3)]), ...
            mkpp(pp.breaks, [6*coef(:,1), 2*coef(:,2)]) ];
    end
end
    
end
