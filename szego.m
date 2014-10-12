classdef szego < cmtobject
% SZEGO class represents a Szego kernel.
%
% S = szego(curve, a)
%   curve - closedcurve object
%   a - point such that f(a) = 0 where f is map to disk.
%
% szego(curve, a, 'name', value, ...)
%   Uses preferences specified by name/value pairs. See szset for valid
%   name/value pairs.
%
% s = theta(S, t) - boundary correspondence function gives angle on the unit
%     circle of the image of a point on curve given by parameter t,
%          th(t) = angle(-1i*phi(t).^2.*tangent(curve, t)),
%     normalized so that th(0) = 0.
%
% sp = thetap(S, t) - derivative of th(t),
%          thp(t) = 2*pi/S(a,a)*abs(phi(t).^2)
%
% t = invtheta(S, s, tol) - inverse of theta via a Newton method with automatic
%     initial guess. Input s is normalized to [0,2*pi), and default tolerance
%     tol is 1e-12.
%
% v = phi(S, t) - Szego kernel interpolation
%          phi(t) := sqrt(abs(tangent(curve, t))).*S(point(curve, t), a)
%     where S(z, a) is the computed Szego kernel.
%
% Trummer, M. "An efficient implementation of a conformal mapping method based
% on the Szego kernel." SIAM J. Numer. Anal., 23(4), 1986.
%
% See also szset.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

properties
    curve             % Closed curve where kernel is defined.
    confCenter = 0    % Point that goes to origin under Riemann map.
    numCollPts        % Number of collocation points.

    phiColl           % Collcation phi values.
    theta0            % Rotation constant.
    Saa               % Szego kernel value S(a,a).
    relResid          % Collocation solution relative residual.

    zPts              % Stored collocation points.
    zTan              % Stored collocation tangents.
    zUnitTan          % Stored collocation unit tangents.
    dtColl = 0        % Collocation differential.

    newtTol           % Newton iteration tolerance.
    beNoisy = false    % Logical; informational output?
end

methods
    function S = szego(curve, varargin)
        % Constructor.
        
        % Make sure defaults are set.
        opts = get(S, szset);
        
        if ~nargin
            return
        end
        
        if ~isa(curve, 'closedcurve')
            error('CMT:InvalidArgument', 'Expected a closedcurve object.')
        end
        
        if nargin > 2
            opts = set(opts, varargin{:});
        end
        
        a = opts.confCenter;
        kerndat = szego.compute_kernel(curve, a, opts);
        knames = fieldnames(kerndat);
        for k = 1:numel(knames)
            S.(knames{k}) = kerndat.(knames{k});
        end

        S.curve = curve;
        S.confCenter = a;
        S.numCollPts = opts.numCollPts;
        S.theta0 = angle(-1i*phi(S, 0)^2*tangent(S.curve, 0));
        S.Saa = sum(abs(S.phiColl.^2))*S.dtColl;
        S.newtTol = opts.newtonTol;
        S.beNoisy = opts.trace;
    end

    function disp(S)
        fprintf('Szego kernel object:\n\n')
        astr = strtrim(evalc('disp(S.confCenter)'));
        fprintf('with a = %s,\n', astr)
        resstr = strtrim(evalc('disp(S.relResid)'));
        fprintf('and kernel solution relative residual\n\t%s,\n', resstr)
        fprintf('computed over curve\n')
        disp(S.curve)
    end

    function A = kerz_stein(S, t)
        % Calculate Kerzmann-Stein kernel.
        t = t(:);
        w = point(S.curve, t);
        wt = tangent(S.curve, t);
        wT = wt./abs(wt);

        z = S.zPts;
        zt = S.zTan;
        zT = S.zUnitTan;
        separation = 10*eps(max(abs(z)));

        function A = KS_by_idx(wi, zi)
            z_w = z(zi).' - w(wi);
            A = sqrt(abs(wt(wi).*zt(zi).'))/(2i*pi).* ...
                (conj(wT(wi)./z_w) - zT(zi).'./z_w);
            A(abs(z_w) < separation) = 0;
        end
        % Function bsxfun will call KS_by_idx with single wi and a vector of zi.
        % Column vector w and row vector z determines shape of resulting
        % matrix A.
        A = bsxfun(@KS_by_idx, (1:numel(w))', 1:S.numCollPts);
    end

    function v = phi(S, t)
        % Calculate scaled Szego kernel.
        v = psi(S, t) - kerz_stein(S, t)*S.phiColl*S.dtColl;
    end

    function y = psi(S, t)
        % Calculate integral equation RHS.
        wt = tangent(S.curve, t(:));
        y = 1i/(2*pi)./sqrt(abs(wt)) .* ...
            conj(wt./(point(S.curve, t(:)) - S.confCenter));
    end

    function th = theta(S, t)
        % Give unit circle correspondence angle.
        th = angle(-1i.*phi(S, t).^2.*tangent(S.curve, t(:))) - S.theta0;
        % KLUDGE: Roundoff allows theta(S, 0) ~= 0. "Fix" this.
        th(t == 0) = 0;
    end

    function [t, output] = invtheta(S, s, tol)
        % Calculate inverse of theta by bisection/Newton method.

        % This should really be more modularised. -- EK

        if nargin < 3 || isempty(tol)
            ntol = S.newtTol;
        else
            ntol = tol;
        end
        trace_label = 'CMT:szego:invtheta';

        % Should check that 1) s is monotonically increasing, and 2) s(end) <
        % 2*pi.
        if any(diff(s) < 0)
            error('CMT:InvalidArgument', 'Input s must be strictly monotonic.')
        end
        if any(s == 2*pi)
            warning('CMT:BadIdea', ...
                ['Passing 2*pi as a value of s may have unpredictible' ...
                ' results. Use 0 instead.'])
        end

        f = @(t, s) s - mod(theta(S, t), 2*pi);
        % Start with a proportionally spaced guess. Hey, you gotta start
        % somewhere, and sometimes you get lucky.
        t = s/(2*pi);

        % Bisection generates initial guess for Newton. We take advantage of
        % the fact that, as a function of t,
        %       f(t, s) := s - theta(S, t)
        % on the proper branch is monotonically decreasing.
        btol = 1e-3;        % Solves problem to order 1e-4 for Newton.
        bmaxiter = 20;      % No runaways!

        % Using convergence rate of bisection method to get an initial partition of
        % [0, length(curve)) so bisection finishes closeish to 5 steps. There ends up
        % being an evaluation at a large number of points, but it's only once.
        nb = max(ceil(1/(2^4*btol)), numel(t));
        if nb > numel(t)
            tt = (0:nb-1)'/nb;
        else
            tt = t;
        end
        th = mod(theta(S, tt), 2*pi);

        % Sign change indicates where zero lives.
        [chg,colk] = find(diff(sign(bsxfun(@minus, s', th))) == -2);
        left = zeros(size(t));
        left(colk) = tt(chg);
        right = ones(size(t));
        right(colk) = tt(chg+1);

        % Bisect.
        done = abs(f(t, s)) < btol;
        biter = 0;
        if S.beNoisy
            fprintf('%s: Starting bisection ...\n', trace_label)
        end
        while ~all(done) && biter < bmaxiter
            biter = biter + 1;

            t(~done) = 0.5*(left(~done) + right(~done));
            fk = f(t(~done), s(~done));
            isneg = fk < 0;
            left(~done) = isneg.*left(~done) + ~isneg.*t(~done);
            right(~done) = isneg.*t(~done) + ~isneg.*right(~done);
            done(~done) = abs(fk) < btol;
        end
        if S.beNoisy
            fprintf('%s: Bisection finished in %d steps.\n', ...
                trace_label, biter)
        end

        % Apply Newton's method.
        nmaxiter = 20;      % Should probably be set somewhere else.

        fval = f(t, s);
        done = abs(fval) < ntol;
        update = double(~done);
        prev_update = nan(size(update));

        niter = 0;
        if S.beNoisy
            fprintf('%s: Starting Newton iteration ...\n', trace_label)
        end
        while ~all(done) && niter < nmaxiter
            niter = niter + 1;

            update(~done) = fval(~done)./thetap(S, t(~done));
            t(~done) = t(~done) + update(~done);

            if all(abs(abs(prev_update(~done)) - abs(update(~done))) <= 100*eps)
                break
            end
            prev_update = update;

            fval(~done) = f(t(~done), s(~done));
            done(~done) = abs(fval(~done)) < ntol;
            update(done) = 0;
        end
        if S.beNoisy
            fprintf('%s: Newton finished in %d steps.\n', ...
                trace_label, niter)
        end

        maxerr = max(abs(fval));
        if S.beNoisy
            fprintf('%s: %d/%d points with |f| > eps, max|f| = %.4e.\n\n', ...
                trace_label, sum(~done), numel(t), max(abs(fval)))
        end
        if maxerr > 100*ntol
            warning('CMT:OutOfTolerance', ...
                ['Newton iteration finished with maxerr=%.6e (>100*tol).\n' ...
                'Check output (maybe szego needs a larger N?).'], ...
                maxerr)
        end

        if nargout > 1
            output.biter = biter;
            output.niter = niter;
            output.maxerr = max(abs(fval));
            output.done = sum(done);
        end
    end

    function thp = thetap(S, t)
        % Derivative of theta.
        thp = 2*pi/S.Saa*abs(phi(S, t(:)).^2);
    end
end

methods(Access=protected, Static)
    function out = compute_kernel(curve, a, opts)
        noisy = opts.trace;
        N = opts.numCollPts;

        % Parameter length should be 1, if not this should be length(curve)/N.
        dt = 1/N;
        t = (0:N-1)'*dt;
        z = point(curve, t);
        zt = tangent(curve, t);
        zT = zt./abs(zt);

        if noisy
            fprintf('\nSzego constructor:\n\n')

            fprintf('  Creating Kerzmann-Stein kernel for %d collocation points...', ...
                N);
        end

        % IpA = I + A, where A is the Kerzmann-Stein kernel.
        IpA = ones(N);
        for j = 2:N
            cols = 1:j-1;
            zc_zj = z(cols) - z(j);
            IpA(j,cols) = (conj(zT(j)./zc_zj) ...
                - zT(cols)./zc_zj) ...
                .* sqrt(abs(zt(j)*zt(cols))) * (dt/(2i*pi));
            IpA(cols,j) = -IpA(j,cols)';
        end

        if noisy
            fprintf('\n')
        end

        y = 1i*sqrt(abs(zt))/(2*pi) .* conj(zT./(z - a));

        if noisy
            fprintf('  Solving for Szego kernel at collocation points using ')
        end

        if strcmp(opts.kernSolMethod, 'auto')
            if N < 2048
                method = 'bs';
            else
                method = 'or';
            end
        else
            method = opts.kernSolMethod;
        end
        if any(strcmp(method, {'backslash', 'bs'}))
            if noisy
                fprintf('backslash...')
            end
            x = IpA\y;
        elseif any(strcmp(method, {'orth_resid', 'or'}))
            if noisy
                fprintf('orthogonal residuals...\n')
            end
            % Need initial guess; get it via interpolation.
            sargs = varargs(opts);
            tmp = szego(curve, sargs{:}, 'numCollPts', 256);
            x = phi(tmp, t);

            % Pass guess to solver.
            x = szego.orthog_resid(IpA, x, y, noisy);
        else
            error('CMT:InvalidArgument', 'Invalid solution method requested.')
        end
        if noisy
            fprintf('\n')
        end

        relresid = norm(y - IpA*x)/norm(y);
        if relresid > 100*eps
            warning('CMT:OutOfTolerance', ...
                'Relative residual %.4g is larger than 100*eps.', ...
                relresid)
        end

        if noisy
            fprintf('  Kernel solution relative residual is %.4e.\n\n', relresid)
        end

        out.phiColl = x;
        out.dtColl = dt;
        out.zPts = z;
        out.zTan = zt;
        out.zUnitTan = zT;
        out.relResid = relresid;
    end

    function x = orthog_resid(A, x, y, noisy)
        % Method of orthogonal residuals. See p. 857 in Trummer.
        if nargin < 4 || isempty(noisy)
            noisy = false;
        end

        if noisy
            fprintf('\n  iter  orth-resid\n')
            fprintf('  ----  ----------\n')
        end

        maxiter = 20;
        omeg = 1;
        xp = x;
        xpp = zeros(size(x));
        yip = y'*y;
        for iter = 1:maxiter
            r = y - A*x;

            rip = r'*r;
            if noisy
                fprintf('   %2d   %.4e\n', iter, sqrt(rip/yip));
            end
            if sqrt(rip/yip) <= eps
                break
            end
            if iter > 1
                omeg = 1/(1 + rip/ripp/omeg);
            end
            ripp = rip;

            x = xpp + omeg*(r + xp - xpp);
            xpp = xp;
            xp = x;
        end
    end
end
    
end
