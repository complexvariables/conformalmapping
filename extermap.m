classdef extermap < scmap
%EXTERMAP Schwarz-Christoffel exterior map object.
%   EXTERMAP(P) constructs a Schwarz-Christoffel exterior map object for
%   the polygon P. The parameter problem is solved using default options
%   for the prevertices and the multiplicative constant.
%
%   EXTERMAP(P,OPTIONS) uses an options structure of the type created by
%   SCMAPOPT in solving the parameter problem.
%
%   EXTERMAP(P,Z) creates a extermap object having the given prevertices
%   Z (the mulitiplicative constant is found automatically).
%   EXTERMAP(P,Z,C) also uses the given constant. An OPTIONS argument
%   can be added, although only the error tolerance will be used.
%
%   EXTERMAP(M), where M is a extermap object, just returns M.
%
%   EXTERMAP(M,P) returns a new extermap object for the polygon P using
%   the options in extermap M. The prevertices of M will be used as the
%   starting guess for the parameter problem of the new map. Thus P
%   should properly be a perturbation of the polygon for M. An OPTIONS
%   structure may also be given to override options in M.
%
%   EXTERMAP(Z,ALPHA) creates a map using the given prevertices and the
%   interior polygon angles described by ALPHA (see POLYGON help). The
%   image polygon is deduced by computing S-C integrals assuming a
%   multiplicative constant of 1. EXTERMAP(Z,ALPHA,C) uses the given
%   constant instead. Note that not every pairing of prevertices and
%   angles produces a single-valued map; you must have SUM((ALPHA-1)./Z)
%   equal to zero. Also, Z is given counterclockwise around the unit
%   circle, but ALPHA should be clockwise with respect to the interior
%   of the polygon.
%
%   See also SCMAPOPT, classes POLYGON, SCMAP.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

%   Copyright 1998-2001 by Toby Driscoll.
%   $Id: extermap.m 129 2001-05-07 15:04:13Z driscoll $

properties
    prevertex
    constant
    qdata
    accuracyVal
end

methods
    function map = extermap(varargin)
        % Initialize with empties
        poly = [];
        alpha = [];
        z = [];
        c = [];
        opt = [];
        qdata = [];
        
        import sctool.*
        
        % Branch based on class of first argument
        switch class(varargin{1})
            case 'extermap'
                oldmap = varargin{1};
                % Continuation of given map to given polygon
                poly = varargin{2};
                opt = scmapopt(oldmap);
                z0 = oldmap.prevertex;
                if length(z0) ~= length(poly)
                    msg = 'Polygon %s must have the same length as that in %s.';
                    error(msg,inputname(2),inputname(1))
                end
                if nargin > 2
                    opt = scmapopt(opt,varargin{3});
                end
                opt = scmapopt(opt,'initial',z0);
                
            case 'polygon'
                poly = varargin{1};
                % Parse optional arguments
                for j = 2:length(varargin)
                    arg = varargin{j};
                    % Each arg is an options struct, z, or c
                    if isa(arg,'struct')
                        opt = arg;
                    elseif length(arg) == length(poly)
                        z = arg;
                        % We will have to flip vertices to get correct orientation
                        z = flipud(z(:));
                    elseif length(arg) == 1
                        c = arg;
                    else
                        msg = 'Unable to parse argument ''%s''.';
                        error(msg,inputname(j+1))
                    end
                end
                
            case 'double'
                % Args are the prevertex vector, then angle vector
                z = varargin{1}(:);
                alpha = varargin{2};
                poly = polygon(NaN*alpha*1i,alpha);
                c = 1;
                % Check residue of integrand to see if compatible
                if abs(sum((alpha-1)./z)) > 1e-8
                    error('Map is not single-valued')
                end
                for j = 3:length(varargin)
                    if isa(varargin{j},'struct')
                        opt = varargin{j};
                    elseif length(varargin{j})==1
                        c = varargin{j};
                    else
                        msg = 'Unable to parse argument ''%s''.';
                        error(msg,inputname(j+1))
                    end
                end
                
            otherwise
                msg = 'Expected ''%s'' to be of class polygon or extermap.';
                error(msg,inputname(1))
        end % switch
        
        
        % Retrieve options
        opt = scmapopt(opt);
        
        % Take actions based on what needs to be filled in
        
        if isempty(z)
            % Find prevertices
            % Apply SCFIX to enforce solver rules
            w = flipud(vertex(poly));
            beta = 1 - flipud(angle(poly));
            [w,beta] = scfix('de',w,beta);
            poly = polygon(flipud(w),1-flipud(beta));
            
            z0 = opt.InitialGuess;
            tol = opt.Tolerance;
            [z,c,qdata] = extermap.deparam(w,beta,z0,opt);
        end
        
        if isempty(qdata)
            % Base quadrature accuracy on given options
            nqpts = ceil(-log10(opt.Tolerance));
            beta = 1 - flipud(angle(poly));
            qdata = scqdata(beta,nqpts);
        end
        
        if isempty(c)
            % Find constant
            w = flipud(vertex(poly));
            beta = 1 - flipud(angle(poly));
            mid = z(1)*exp(i*angle(z(2)/z(1))/2);
            I = extermap.dequad(z(1),mid,1,z,beta,qdata) - ...
                extermap.dequad(z(2),mid,2,z,beta,qdata);
            c = diff(w(1:2))/I;
        end
        
        map = map@scmap(poly,opt);
        
        map.prevertex = z;
        map.constant = c;
        map.qdata = qdata;
        
        % If the polygon was not known, find it from the map
        if any(isnan(vertex(poly)))
            poly = forwardpoly(map);
            map.scmap = scmap(poly,opt);
        end
        
        % Now fill in apparent accuracy
        map.accuracyVal = accuracy(map);
    end
    
    function acc = accuracy(M)
        %ACCURACY Apparent accuracy of Schwarz-Christoffel disk exterior map.
        %   ACCURACY(M) estimates the accuracy of the Schwarz-Christoffel
        %   exterior map M. The technique used is to compare the differences
        %   between successive vertices to the integral between the
        %   corresponding prevertices, and return the maximum.
        %
        %   See also EXTERMAP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: accuracy.m 7 1998-05-10 04:37:19Z tad $
        
        % If an accuracy has been assigned, don't question it
        if ~isempty(M.accuracyVal)
            acc = M.accuracyVal;
            return
        end
        
        % Get data for low-level functions
        p = M.polygon;
        w = flipud(vertex(p));
        beta = flipud(1 - angle(p));
        z = M.prevertex;
        c = M.constant;
        qdata = M.qdata;
        n = length(w);
        
        % Test accuracy by integrating between consecutive finite prevertices, and
        % comparing to differences of vertices.
        
        idx = find(~isinf(w));
        
        % Two columns hold endpoint indices for integrations
        idx = [idx(1:end) idx([2:end 1])];
        
        % Find midpoints that are halfway between in angular sense
        dtheta = mod(angle(z(idx(:,2))./z(idx(:,1))),2*pi);
        mid = z(idx(:,1)).*exp(i*dtheta/2);
        
        % Do the integrations
        I = extermap.dequad(z(idx(:,1)),mid,idx(:,1),z,beta,qdata) - ...
            extermap.dequad(z(idx(:,2)),mid,idx(:,2),z,beta,qdata);
        
        acc = max(abs( c*I - diff(w([1:end 1])) ));
    end
    
    function g = capacity(map)
        g = abs(map.constant);
    end
    
    function out = char(f)
        %CHAR   Pretty-print a Schwarz-Christoffel exterior map.
        
        %   Copyright 1998-2001 by Toby Driscoll.
        %   $Id: char.m 161 2001-07-20 14:32:59Z driscoll $
        
        p = f.polygon;
        w = vertex(p);
        alpha = angle(p);
        n = length(w);
        z = flipud(f.prevertex);
        c = f.constant;
        
        L = cell(2,1);
        L{1}='      vertex              alpha         prevertex             arg/pi';
        L{2}=' -----------------------------------------------------------------------';
        u = real(w);
        v = imag(w);
        x = real(z);
        y = imag(z);
        ang = angle(z)/pi;
        ang(ang<=0) = ang(ang<=0) + 2;
        fmt = ' %8.5f %c %7.5fi    %8.5f   %8.5f %c %7.5fi    %14.12f';
        for j = 1:length(w)
            if v(j) < 0
                s1 = '-';
            else
                s1 = '+';
            end
            if y(j) < 0
                s2 = '-';
            else
                s2 = '+';
            end
            L{end+1}=sprintf(fmt,u(j),s1,abs(v(j)),alpha(j),x(j),s2,abs(y(j)),ang(j));
        end
        
        L{end+1} = ' ';
        if imag(c) < 0
            s = '-';
        else
            s = '+';
        end
        L{end+1} = sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c)));
        L{end+1} = sprintf('  Logarithmic capacity = %.8g',abs(c));
        L{end+1} = sprintf('  Apparent accuracy is %.2e',f.accuracyVal);
        L{end+1} = ' ';
        
        out = L;
    end
    
    function out = display(M)
        %DISPLAY Display parameters of a Schwarz-Christoffel exterior map.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: display.m 118 2001-05-03 21:18:27Z driscoll $
        
        p = M.polygon;
        w = vertex(p);
        alpha = angle(p);
        n = length(w);
        z = flipud(M.prevertex);
        c = M.constant;
        
        L = { ' '; '  extermap object:'; ' ' };
        L{4}='      vertex              alpha         prevertex             arg/pi';
        L{5}=' -----------------------------------------------------------------------';
        u = real(w);
        v = imag(w);
        x = real(z);
        y = imag(z);
        ang = angle(z)/pi;
        ang(ang<=0) = ang(ang<=0) + 2;
        fmt = ' %8.5f %c %7.5fi    %8.5f   %8.5f %c %7.5fi    %14.12f';
        for j = 1:length(w)
            if v(j) < 0
                s1 = '-';
            else
                s1 = '+';
            end
            if y(j) < 0
                s2 = '-';
            else
                s2 = '+';
            end
            L{end+1}=sprintf(fmt,u(j),s1,abs(v(j)),alpha(j),x(j),s2,abs(y(j)),ang(j));
        end
        
        L{end+1} = ' ';
        if imag(c) < 0
            s = '-';
        else
            s = '+';
        end
        L{end+1} = sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c)));
        L{end+1} = sprintf('  Logarithmic capacity = %.8g',abs(c));
        L{end+1} = ' ';
        L{end+1} = sprintf('  Apparent accuracy is %.2e',M.accuracyVal);
        L{end+1} = ' ';
        
        
        if nargout==0
            fprintf('%s\n',L{:})
        else
            out = L;
        end
    end
    
    function wp = eval(M,zp,tol)
        %EVAL Evaluate Schwarz-Christoffel exterior map at points.
        %   EVAL(M,ZP) evaluates the Schwarz-Christoffel map M at the points
        %   ZP in the unit disk. The default tolerance of M is used.
        %
        %   EVAL(M,ZP,TOL) attempts to give an answer accurate to TOL. If TOL
        %   is less than the accuracy of M, this is unlikely to be met.
        %
        %   See also EXTERMAP, EVALINV.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: eval.m 7 1998-05-10 04:37:19Z tad $
        
        if nargin < 3
            qdata = M.qdata;
        else
            qdata = tol;
        end
        
        p = M.polygon;
        w = flipud(vertex(p));
        beta = flipud(1 - angle(p));
        n = length(w);
        
        wp = NaN*zp;
        idx = abs(zp) <= 1+eps;
        wp(idx) = M.demap(zp(idx),w,beta,M.prevertex,M.constant,qdata);
    end
    
    function fp = evaldiff(M,zp)
        %EVALDIFF Derivative of Schwarz-Christoffel exterior map at points.
        %   EVALDIFF(M,ZP) computes the derivative of the Schwarz-Christoffel
        %   exterior map M at the points ZP.
        %
        %   See also EXTERMAP, EVAL.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: evaldiff.m 7 1998-05-10 04:37:19Z tad $
        
        z = M.prevertex;
        c = M.constant;
        beta = flipud(1 - angle(M.polygon));
        
        fp = M.dederiv(zp,z,beta,c);
    end
    
    function zp = evalinv(M,wp,tol,z0)
        %EVALINV Invert Schwarz-Christoffel exterior map at points.
        %   EVALINV(M,WP) evaluates the inverse of the Schwarz--Christoffel map M
        %   at the points WP in the polygon. The default tolerance of M is used.
        %
        %   EVALINV(M,WP,TOL) attempts to give an answer accurate to TOL. If TOL
        %   is smaller than the accuracy of M, this is unlikely to be met.
        %
        %   EVALINV(M,WP,TOL,Z0) uses given starting points. Z0 must be either
        %   the same size as WP or a complex scalar (to be expanded to that
        %   size). It is used for the starting approximation to the inverse
        %   image of WP. The starting guess need not be close to the correct
        %   answer; however, the straight line segment between WP(K) and the
        %   forward image of Z0(K) must lie entirely inside the polygon, for
        %   each K.
        %
        %   See also EXTERMAP, EVAL.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: evalinv.m 7 1998-05-10 04:37:19Z tad $
        
        % Assign empties to missing input args
        if nargin < 4
            z0 = [];
            if nargin < 3
                tol = [];
            end
        end
        
        % Check inputs/supply defaults
        if ~isempty(tol)
            % Either a scalar tolerance or a qdata matrix was passed
            qdata = tol;
            if length(tol) > 1
                tol = 10^(-size(qdata,1));
            end
        else
            qdata = M.qdata;
            tol = M.accuracyVal;
        end
        
        if ~isempty(z0)
            if length(z0) == 1
                z0 = repmat(z0,size(wp));
            elseif any(size(z0) ~= size(wp))
                msg = 'Argument %s must be a complex scalar or the same size as %s.';
                error(sprintf(msg,inputname(z0),inputname(1)));
            end
        end
        
        p = M.polygon;
        w = flipud(vertex(p));
        beta = flipud(1 - angle(p));
        n = length(w);
        
        zp = NaN*wp;
        %idx = logical(~isinpoly(wp,p));
        idx = logical(ones(size(wp)));
        zp(idx) = M.deinvmap(wp(idx),w,beta,M.prevertex,M.constant,qdata,z0,[0 tol]);
    end
    
    function varargout = feval(varargin)
        %FEVAL   Equivalent to EVAL.
        
        if nargout
            varargout = cell(1,nargout);
            [varargout{:}] = eval(varargin{:});
        else
            varargout{1} = eval(varargin{:});
        end
    end
    
    function p = forwardpoly(map)
        %   Given an extermap M, FORWARDPOLY(M) returns the polygon that is
        %   formed using the prevertices, angles, and quadrature data of that
        %   map. If the prevertices were found from the solution of a
        %   parameter problem, then the result should agree closely with the
        %   original polygon that was supplied.
        
        %   Copyright (c) 1998 by Toby Driscoll.
        %   $Id: forwardpoly.m 39 1998-07-01 17:40:15Z tad $
        
        z = map.prevertex;
        alpha = flipud(angle(map.polygon));
        c = map.constant;
        
        n = length(z);
        
        % Since there is no parameter problem, use high accuracy in quadrature.
        qdata = scqdata(1-alpha,16);
        
        % Midpoints of integration
        theta = rem(angle(z(n)) + angle(z/z(n))+2*pi,2*pi);
        theta(end) = 2*pi;
        mid = exp(i*(theta(1:n-1)+theta(2:n))/2);
        
        % Integrations
        I = map.dequad(z(1:n-1),mid,1:n-1,z,1-alpha,qdata) - ...
            map.dequad(z(2:n),mid,2:n,z,1-alpha,qdata);
        
        % Deduce vertices
        w = c*cumsum([0;I]);
        
        p = w.polygon;
    end
    
    function varargout = get(map,varargin)
        %GET    Get map parameters.
        %   [VAL1,VAL2,...] = GET(F,'PROP1','PROP2',...) returns the values of the
        %   map F corresponding to the requested properties. Valid properties
        %   are:
        %
        %       polygon, options, prevertex, constant
        
        % Copyright 1999-2003 by Toby Driscoll.
        % $Id: get.m 236 2003-01-15 15:29:14Z driscoll $
        
        for j = 1:length(varargin)
            switch lower(varargin{j}(1:min(3,length(varargin{j}))))
                case 'pol'
                    varargout{j} = map.scmap.polygon;
                case 'opt'
                    varargout{j} = map.scmap.options;
                case 'pre'
                    varargout{j} = map.prevertex;
                case 'con'
                    varargout{j} = map.constant;
                otherwise
                    warning(sprintf('Property ''%s'' not recognized.\n',varargin{j}))
                    varargout{j} = [];
            end
        end
    end
    
    function M = mtimes(M,c)
        %   Scale the image of a map by a complex constant.
        
        %   Copyright (c) 1998 by Toby Driscoll.
        %   $Id: mtimes.m 33 1998-06-29 22:35:40Z tad $
        
        % May need to swap arguments
        if isa(M,'double') & isa(c,'extermap')
            tmp = M;
            M = c;
            c = tmp;
        end
        
        M.constant = c*M.constant;
        M.scmap = c*M.scmap;
    end
    
    function v = parameters(M)
        %PARAMETERS Return a structure of the Schwarz-Christoffel map parameters.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: parameters.m 7 1998-05-10 04:37:19Z tad $
        
        v.prevertex = flipud(M.prevertex);
        v.constant = M.constant;
    end
    
    function [h,r,theta] = plot(M,varargin)
        %PLOT Visualize a Schwarz-Christoffel exterior map.
        %   PLOT(M) plots the polygon associated with the Schwarz-Christoffel
        %   exterior map M and the images of ten evenly spaced circles and radii
        %   under the S-C transformation.
        %
        %   PLOT(M,NR,NTHETA) plots the images of NR circles and NTHETA radii.
        %
        %   PLOT(M,R,THETA) plots the circles at radii given by the entries of R
        %   and radii at the angles specified in THETA.
        %
        %   PLOT(M,TOL) or PLOT(M,NR,NTHETA,TOL) or PLOT(M,R,THETA,TOL)
        %   computes the map with accuracy roughly TOL. Normally TOL defaults to
        %   1e-4 or the accuracy of M, whichever is greater.
        %
        %   See also EXTERMAP, MAP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: plot.m 7 1998-05-10 04:37:19Z tad $
        
        p = M.polygon;
        w = flipud(vertex(p));
        beta = flipud(1 - angle(p));
        z = M.prevertex;
        c = M.constant;
        n = length(w);
        
        if nargin == 1
            [a1,a2,a3] = M.deplot(w,beta,z,c);
        elseif length(varargin) == 1
            % Tolerance given only
            [a1,a2,a3] = M.deplot(w,beta,z,c,10,10,ceil(-log10(varargin{1})));
        elseif length(varargin) == 2
            % R, theta given only
            [a1,a2,a3] = M.deplot(w,beta,z,c,varargin{1},varargin{2});
        else
            % All given
            nqpts = ceil(-log10(varargin{3}));
            [a1,a2,a3] = M.deplot(w,beta,z,c,varargin{1},varargin{2},nqpts);
        end
        
        if nargout > 0
            h = a1;
            r = a2;
            theta = a3;
        end
    end
end

methods(Hidden,Static)
    function fprime = dederiv(zp,z,beta,c)
        %DEDERIV Derivative of the exterior map.
        %   DEDERIV(ZP,Z,BETA,C) returns the derivative at the points of ZP of the
        %   Schwarz-Christoffel exterior map defined by Z, BETA, and C.
        %
        %   See also DEPARAM, DEMAP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: dederiv.m 7 1998-05-10 04:37:19Z tad $
        
        % Support old syntax
        if nargin < 4
            c = 1;
        end
        
        z = z(:);
        beta = [beta(:);-2];
        zprow = zp(:).';
        fprime = zeros(size(zp));
        
        npts = length(zp(:));
        terms = 1 - zprow(ones(length(z),1),:)./z(:,ones(npts,1));
        terms(length(z)+1,:) = zprow;
        fprime(:) = c*exp(sum(log(terms).*beta(:,ones(npts,1))));
    end
    
    function dedisp(w,beta,z,c)
        %DEDISP Display results of Schwarz-Christoffel exterior parameter problem.
        %   DEDISP(W,BETA,Z,C) displays the results of DEPARAM in a pleasant
        %   way.
        %
        %   See also DEPARAM, DEPLOT.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: dedisp.m 7 1998-05-10 04:37:19Z tad $
        
        disp(' ')
        disp('          w               beta              z               arg(z)/pi')
        disp(' -----------------------------------------------------------------------')
        u = real(w);
        v = imag(w);
        x = real(z);
        y = imag(z);
        ang = angle(z)/pi;
        ang(ang<=0) = ang(ang<=0) + 2;
        for j = 1:length(w)
            if v(j) < 0
                s1 = '-';
            else
                s1 = '+';
            end
            if y(j) < 0
                s2 = '-';
            else
                s2 = '+';
            end
            disp(sprintf(' %8.5f %c %7.5fi    %8.5f   %8.5f %c %7.5fi    %14.12f',...
                u(j),s1,abs(v(j)),beta(j),x(j),s2,abs(y(j)),ang(j)));
            
        end
        disp(' ')
        if imag(c) < 0
            s = '-';
        else
            s = '+';
        end
        disp(sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c))))
    end
    
    function zp = deinvmap(wp,w,beta,z,c,qdat,z0,options)
        %DEINVMAP Schwarz-Christoffel exterior inverse map.
        %   DEINVMAP(WP,W,BETA,Z,C,TOL) computes the inverse of the
        %   Schwarz-Christoffel exterior map (i.e., from the exterior of a
        %   polygon to the disk) at the points given in vector WP.  The other
        %   arguments are as in DEPARAM.  TOL is a scalar tolerance, or a
        %   quadrature-data matrix QDAT as returned by SCQDATA, or may be
        %   omitted.
        %
        %   The default algorithm is to solve an ODE in order to obtain a fair
        %   approximation for ZP, and then improve ZP with Newton iterations.
        %   The ODE solution at WP requires a vector Z0 whose forward image W0
        %   is such that for each j, the line segment connecting WP(j) and W0(j)
        %   lies inside the polygon.  By default Z0 is chosen by a fairly robust
        %   automatic process.  Using a parameter (see below), you can choose to
        %   use either an ODE solution or Newton iterations exclusively.
        %
        %   DEINVMAP(WP,W,BETA,Z,C,TOL,Z0) has two interpretations.  If the ODE
        %   solution is being used, Z0 overrides the automatic selection of
        %   initial points.  (This can be handy in convex polygons, where the
        %   choice of Z0 is trivial.)  Otherwise, Z0 is taken as an initial
        %   guess to ZP.  In either case, if length(Z0)==1, the value Z0 is used
        %   for all elements of WP; otherwise, length(Z0) should equal
        %   length(WP).
        %
        %   DEINVMAP(WP,W,BETA,Z,C,TOL,Z0,OPTIONS) uses a vector of parameters
        %   that control the algorithm.  See SCINVOPT.
        %
        %   See also SCINVOPT, DEPARAM, DEMAP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: deinvmap.m 279 2007-05-14 20:14:24Z driscoll $
        
        n = length(w);
        z = z(:);
        beta = beta(:);
        zp = zeros(size(wp));
        wp = wp(:);
        lenwp = length(wp);
        import sctool.*
        
        if nargin < 8
            options = [];
            if nargin < 7
                z0 = [];
                if nargin < 6
                    qdat = [];
                end
            end
        end
        
        [ode,newton,tol,maxiter] = scinvopt(options);
        
        if isempty(qdat)
            qdat = tol;
        end
        
        if length(qdat)==1
            qdat = scqdata(beta,max(ceil(-log10(qdat)),2));
        end
        
        done = zeros(size(wp));
        % First, trap all points indistinguishable from vertices, or they will cause
        % trouble.
        % Modified 05/14/2007 to work around bug in matlab 2007a.
        for j=1:n
            idx = find(abs(wp-w(j)) < 3*eps);
            zp(idx) = z(j);
            done(idx) = 1;
        end
        lenwp = lenwp - sum(done);
        if lenwp==0, return, end
        
        % ODE
        if ode
            if isempty(z0)
                % Pick a value z0 (not a singularity) and compute the map there.
                map = @(zp) extermap.demap(zp,w,beta,z,c,qdat);
                [z0,w0] = sctool.findz0('de',wp(~done),map,w,beta,z,c,qdat);
            else
                w0 = extermap.demap(z0,w,beta,z,c,qdat);
                if length(z0)==1 & lenwp > 1
                    z0 = z0(:,ones(lenwp,1)).';
                    w0 = w0(:,ones(lenwp,1)).';
                end
                w0 = w0(~done);
                z0 = z0(~done);
            end
            
            % Use relaxed ODE tol if improving with Newton.
            odetol = max(tol,1e-3*(newton));
            
            % Rescale dependent coordinate
            scale = (wp(~done) - w0(:));
            
            % Solve ODE
            z0 = [real(z0);imag(z0)];
            odefun = @(w,y) extermap.deimapfun(w,y,scale,z,beta,c);
            [t,y] = ode23(odefun,[0,0.5,1],z0,odeset('abstol',odetol));
            [m,leny] = size(y);
            zp(~done) = y(m,1:lenwp)+sqrt(-1)*y(m,lenwp+1:leny);
            out = abs(zp) > 1;
            zp(out) = sign(zp(out));
        end
        
        % Newton iterations
        if newton
            if ~ode
                zn = z0(:);
                if length(z0)==1 && lenwp > 1
                    zn = zn(:,ones(lenwp,1));
                end
                zn(done) = zp(done);
            else
                zn = zp(:);
            end
            
            wp = wp(:);
            k = 0;
            while ~all(done) & k < maxiter
                F = wp(~done) - extermap.demap(zn(~done),w,beta,z,c,qdat);
                m = length(F);
                dF = c*(zn(~done).').^(-2) .* exp(sum(beta(:,ones(m,1)) .* ...
                    log(1-(zn(~done,ones(n,1)).')./z(:,ones(m,1)))));
                zn(~done) = zn(~done) + F(:)./dF(:);
                done(~done) = (abs(F)< tol);
                k = k+1;
            end
            if any(abs(F)> tol)
                str = sprintf('Check solution; maximum residual = %.3g\n',max(abs(F)));
                warning(str)
            end
            zp(:) = zn;
        end
    end

    function zdot = deimapfun(wp,yp,scale,z,beta,c)
        %   Used by DEINVMAP for solution of an ODE.

        %   Copyright 1998 by Toby Driscoll.
        %   $Id: deimapfun.m 7 1998-05-10 04:37:19Z tad $

        lenyp = length(yp);
        lenzp = lenyp/2;
        zp = yp(1:lenzp)+sqrt(-1)*yp(lenzp+1:lenyp);

        f = scale./extermap.dederiv(zp,z,beta,c);
        zdot = [real(f);imag(f)];
    end

    function wp = demap(zp,w,beta,z,c,qdat)
        %DEMAP  Schwarz-Christoffel exterior map.
        %   DEMAP(ZP,W,BETA,Z,C,QDAT) computes the values of the Schwarz-
        %   Christoffel exterior map at the points in vector ZP. The arguments
        %   W, BETA, Z, C, and QDAT are as in DEPARAM.  DEMAP returns a vector
        %   the same size as ZP.
        %
        %   DEMAP(ZP,W,BETA,Z,C,TOL) uses quadrature data intended to give an
        %   answer accurate to within TOL.
        %
        %   See also DEPARAM, DEPLOT, DEINVMAP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: demap.m 7 1998-05-10 04:37:19Z tad $
        
        if isempty(zp)
            wp = [];
            return
        end
        
        n = length(w);
        beta = beta(:);
        z = z(:);
        
        % Quadrature data and error tolerance
        if nargin < 6
            tol = 1e-8;
            qdat = scqdata(beta,8);
        elseif length(qdat)==1
            tol = qdat;
            qdat = scqdata(beta,max(ceil(-log10(tol)),8));
        else
            tol = 10^(-size(qdat,1));
        end
        
        shape = size(zp);
        zp = zp(:);
        zprow = zp.';
        p = length(zp);
        wp = zeros(p,1);
        
        % For each point in zp, find nearest prevertex.
        [dist,sing] = min(abs(zprow(ones(n,1),:) - z(:,ones(1,p))));
        sing = sing(:);				% indices of prevertices
        
        % Screen out images of prevertices
        vertex = (dist(:) < tol);
        wp(vertex) = w(sing(vertex));
        zero = abs(zp) < tol;
        wp(zero) = Inf;
        vertex = vertex | zero;
        
        % zs = the starting singularities
        zs = z(sing);
        
        wp(~vertex) = w(sing(~vertex));
        
        % Must be careful about the singularity at the origin, since the
        % quadrature routine doesn't pay attention to the right endpoint.
        
        abszp = abs(zp); 			% dist to sing at 0
        
        % Integrate for the rest.
        unf = ~vertex;
        zold = zs(unf);
        while any(unf)
            % How far can the integration go?
            dist = min(1,2*abszp(unf)./abs(zp(unf)-zold));
            % New integration endpoints
            znew = zold + dist.*(zp(unf)-zold);
            wp(unf) = wp(unf) + ...
                c*extermap.dequad(zold,znew,sing(unf),z,beta,qdat);
            
            unf = (dist<1);
            zold = znew(unf);
            sing(:) = 0;
        end
        
        wp = reshape(wp,shape);
    end
    
    function [z,c,qdat] = deparam(w,beta,z0,options)
        %DEPARAM Schwarz-Christoffel exterior parameter problem.
        %   [Z,C,QDAT] = DEPARAM(W,BETA) solves the Schwarz-Christoffel mapping
        %   parameter problem with a disk as fundamental domain and the exterior
        %   of the polygon specified by W as the target. W must be a vector of
        %   the vertices of the polygon, specified in clockwise order, and BETA
        %   should be a vector of the turning angles of the polygon; see
        %   SCANGLES for details. If successful, DEPARAM will return Z, a vector
        %   of the pre-images of W; C, the multiplicative constant of the
        %   conformal map; and QDAT, an optional matrix of quadrature data used
        %   by some of the other S-C routines.
        %
        %   [Z,C,QDAT] = DEPARAM(W,BETA,Z0) uses Z0 as an initial guess for Z.
        %
        %   [Z,C,QDAT] = DEPARAM(W,BETA,TOL) attempts to find an answer within
        %   tolerance TOL. (Also see SCPAROPT.)
        %
        %   [Z,C,QDAT] = DEPARAM(W,BETA,Z0,OPTIONS) uses a vector of control
        %   parameters. See SCPAROPT.
        %
        %   See also SCPAROPT, DRAWPOLY, DEDISP, DEPLOT, DEMAP, DEINVMAP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: deparam.m 129 2001-05-07 15:04:13Z driscoll $
        
        n = length(w); 				% no. of vertices
        w = w(:);
        beta = beta(:);
        import sctool.*
        
        % Set up defaults for missing args
        if nargin < 4
            options = [];
            if nargin < 3
                z0 = [];
            end
        end
        
        err = sccheck('de',w,beta);
        if err==1
            fprintf('Use SCFIX to make polygon obey requirements\n')
            error(' ')
        end
        
        [trace,tol,method] = parseopt(options);
        if length(z0)==1
            tol = z0;
            z0 = [];
        end
        nqpts = max(ceil(-log10(tol)),2);
        qdat = scqdata(beta,nqpts); 		% quadrature data
        
        if n==2					% it's a slit
            z = [-1;1];
            
        else
            % Set up normalized lengths for nonlinear equations
            len = abs(diff(w([n,1:n])));
            nmlen = abs(len(3:n-1)/len(2));
            
            % Set up initial guess
            if isempty(z0)
                y0 = zeros(n-1,1);
            else
                z0 = z0/z0(n);			% Fix z0(n)=1.
                th = angle(z0(:));
                th(th<=0) = th(th<=0) + 2*pi;
                dt = diff([0;th(1:n-1);2*pi]);
                y0 = log(dt(1:n-1)./dt(2:n));
            end
            
            % Solve nonlinear system of equations:
            
            % package data
            fdat = {n,beta,nmlen,qdat};
            % set options
            opt = zeros(16,1);
            opt(1) = trace;
            opt(2) = method;
            opt(6) = 100*(n-3);
            opt(8) = tol;
            opt(9) = min(eps^(2/3),tol/10);
            opt(12) = nqpts;
            try
                [y,termcode] = nesolve(@extermap.depfun,y0,opt,fdat);
            catch
                % Have to delete the "waitbar" figure if interrupted
                close(findobj(allchild(0),'flat','Tag','TMWWaitbar'));
                error(lasterr)
            end
            if termcode~=1
                warning('Nonlinear equations solver did not terminate normally.')
            end
            
            % Convert y values to z
            cs = cumsum(cumprod([1;exp(-y)]));
            theta = 2*pi*cs/cs(n);
            z = ones(n,1);
            z(1:n-1) = [exp(1i*theta(1:n-1))];
        end
        
        % Determine scaling constant
        mid = z(1)*exp(1i*angle(z(2)/z(1))/2);
        c = (w(2) - w(1)) / (extermap.dequad(z(1),mid,1,z,beta,qdat)-...
            extermap.dequad(z(2),mid,2,z,beta,qdat));
    end

    function F = depfun(y,fdat)
        %   Returns residual for solution of nonlinear equations.
        %   $Id: depfun.m 268 2003-05-08 18:08:25Z driscoll $

        [n,beta,nmlen,qdat] = deal(fdat{:});

        % Transform y (unconstr. vars) to z (prevertices)
        cs = cumsum(cumprod([1;exp(-y)]));
        theta = 2*pi*cs(1:n-1)/cs(n);
        z = ones(n,1);
        z(1:n-1) = exp(1i*theta);

        % Compute the integrals
        mid = exp(1i*(theta(1:n-2)+theta(2:n-1))/2);

        % We can use the same quadrature as for the interior map, because the abs
        % value of the integrand on the unit circle is not affected by the z^{-2}
        % term.
        ints = sctool.dabsquad(z(1:n-2),mid,1:n-2,z,beta,qdat) + ...
            sctool.dabsquad(z(2:n-1),mid,2:n-1,z,beta,qdat);

        if any(ints==0)
            % Singularities were too crowded in practice.
            warning('Severe crowding')
        end

        % Compute equation residual values.
        if n > 3
            F = abs(ints(2:n-2))/abs(ints(1)) - nmlen;
        else
            F = [];
        end

        % Compute residue.
        res = -sum(beta./z)/ints(1);

        F = [F;real(res);imag(res)];
    end
    
    function [H,R2,THETA] = deplot(w,beta,z,c,R,theta,options)
        %DEPLOT Image of polar grid under Schwarz-Christoffel exterior map.
        %   DEPLOT(W,BETA,Z,C) will adaptively plot the images under the
        %   Schwarz-Christoffel exterior map of ten evenly spaced circles and
        %   rays in the unit disk.  The arguments are as in DEPARAM.
        %
        %   DEPLOT(W,BETA,Z,C,M,N) will plot images of M evenly spaced circles
        %   and N evenly spaced rays.
        %
        %   DEPLOT(W,BETA,Z,C,R,THETA) will plot images of circles whose radii
        %   are given in R and rays whose arguments are given in THETA.  Either
        %   argument may be empty.
        %
        %   DEPLOT(W,BETA,Z,C,R,THETA,OPTIONS) allows customization of DEPLOT's
        %   behavior.  See SCPLTOPT.
        %
        %   H = DEPLOT(W,BETA,Z,C,...) returns a vector of handles to all the
        %   curves drawn in the interior of the polygon.  [H,R,THETA] =
        %   DEPLOT(W,BETA,Z,C,...) also returns the moduli and arguments of the
        %   curves comprising the grid.
        %
        %   See also SCPLTOPT, DEPARAM, DEMAP, DEDISP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: deplot.m 226 2003-01-08 16:59:17Z driscoll $
        
        beta = beta(:);
        z = z(:);
        import sctool.*
        
        % Parse input
        if nargin < 7
            options = [];
            if nargin < 6
                theta = [];
                if nargin < 5
                    R = [];
                end
            end
        end
        
        % Empty arguments default to 10
        if isempty([R(:);theta(:)])
            R = 10;
            theta = 10;
        end
        
        % Integer arguments must be converted to specific values
        if (length(R)==1) && (R == round(R))
            m = R+2;
            R = fliplr(linspace(.25,1,m));
            R([1,m]) = [];
        end
        if (length(theta)==1) && (theta == round(theta))
            m = theta+1;
            theta = linspace(0,2*pi,m);
            theta(m) = [];
        end
        
        autoscale = strcmp(get(gca,'xlimmode'),'auto') && ...
            strcmp(get(gca,'ylimmode'),'auto');
        autoscale = autoscale | ~ishold;
        
        % Prepare figure
        fig = gcf;
        figure(fig);
        
        % If figure has two tagged axes, draw in both
        ax = [findobj(fig,'tag','PhysicalAxes');...
            findobj(fig,'tag','CanonicalAxes')];
        if length(ax)==2
            draw2 = true;
            vis = get(ax,'vis');
        else
            draw2 = false;
            ax = gca;
            vis = {'on'};
        end
        
        % Prepare axes
        axes(ax(1))
        turn_off_hold = ~ishold;
        hp = plotpoly(w,beta);
        set(hp,'vis',vis{1})
        hold on
        % For now, there is no need to draw the canonical domain. This is an
        % awkward truce with the GUI.
        
        axlim = axis;
        if autoscale
            axlim(1:2) = axlim(1:2) + 0.25*diff(axlim(1:2))*[-1,1];
            axlim(3:4) = axlim(3:4) + 0.25*diff(axlim(3:4))*[-1,1];
            axis(axlim);
        end
        drawnow
        
        % Enlarge clipping axes to avoid cutting off arcs
        axlim(1:2) = axlim(1:2) + 0.25*diff(axlim(1:2))*[-1,1];
        axlim(3:4) = axlim(3:4) + 0.25*diff(axlim(3:4))*[-1,1];
        
        % Drawing parameters
        [nqpts,minlen,maxlen,maxrefn] = scpltopt(options);
        qdat = scqdata(beta,nqpts);
        len = max(diff(get(ax(1),'xlim')),diff(get(ax(1),'ylim')));
        minlen = len*minlen;
        maxlen = len*maxlen;
        
        color = 'k';
        
        % Plot circles...
        linh = gobjects(length(R),2);
        for j = 1:length(R)
            % Start with evenly spaced theta
            tp = linspace(0,2*pi,20)';
            new = true(length(tp),1);
            wp = NaN(size(new));
            
            % The individual points will be shown as they are found
            if verLessThan('matlab','8.4')
                linh(j,1) = line(NaN,NaN,'parent',ax(1),'color',color,'vis',vis{1},...
                    'linestyle','none','marker','.','markersize',7,'erasemode','none');
                if draw2
                    linh(j,2) = line(NaN,NaN,'parent',ax(2),'color',color,'vis',vis{2},...
                        'linestyle','none','marker','.','markersize',7,'erasemode','none');
                end
            else
                
                linh(j,1) = animatedline('parent',ax(1),'color',color,'vis',vis{1},...
                    'linestyle','none','marker','.','markersize',7);
                if draw2
                    linh(j,2) = animatedline('parent',ax(2),'color',color,'vis',vis{2},...
                        'linestyle','none','marker','.','markersize',7);
                end
            end
            
            % Adaptively refine theta to make smooth curve
            iter = 0;
            while (any(new)) && (iter < maxrefn)
                drawnow
                zp = R(j)*exp(1i*tp(new));
                neww = extermap.demap(zp,w,beta,z,c,qdat);
                wp(new) = neww;
                iter = iter + 1;
                
                % Update the points to show progress
                if verLessThan('matlab','8.4')
                    set(linh(j,1),'xdata',real(wp),'ydata',imag(wp))
                    if draw2
                        set(linh(j,2),'xdata',R(j)*cos(tp),'ydata',R(j)*sin(tp))
                    end
                else
                    addpoints(linh(j,1),real(wp(new)),imag(wp(new)))
                    if draw2
                        addpoints(linh(j,1),R(j)*cos(tp(new)),R(j)*sin(tp(new)))
                    end
                end
                drawnow
                
                % Add points to zp where necessary
                [tp,wp,new] = scpadapt(tp,wp,minlen,maxlen,axis);
                
            end
            % Set the lines to be solid
            if verLessThan('matlab','8.4')
                set(linh(j,1),'erasemode','back')
                set(linh(j,1),'marker','none','linestyle','-','user',R(j)*exp(1i*tp))
                if draw2
                    % Replace the points with (hopefully) a smooth circle
                    tp = linspace(0,2*pi,101);
                    set(linh(j,2),'erasemode','back')
                    set(linh(j,2),'marker','none','linestyle','-',...
                        'xdata',R(j)*cos(tp),'ydata',R(j)*sin(tp))
                end
                
            else
                
                clearpoints(linh(j,1))
                addpoints(linh(j,1),real(wp),imag(wp));
                set(linh(j,1),'marker','none','linestyle','-','user',R(j)*exp(1i*tp))
                if draw2
                    % Replace the points with (hopefully) a smooth circle
                    tp = linspace(0,2*pi,361);
                    clearpoints(linh(j,2))
                    addpoints(linh(j,2),R(j)*cos(tp),R(j)*sin(tp))
                    set(linh(j,2),'marker','none','linestyle','-')
                end
            end
            
            drawnow
            
        end
        
        
        % Plot radii...
        linh1 = linh;
        linh = gobjects(length(theta),2);
        for j = 1:length(theta)
            Rp = [0 linspace(.2,1,14)]';
            zp = Rp*exp(1i*theta(j));
            new = true(length(zp),1);
            wp = NaN(size(new));
            new(1) = 0;
            wp(1) = Inf;
            
            % The individual points will be shown as they are found
            if verLessThan('matlab','8.4')
                linh(j,1) = line(NaN,NaN,'parent',ax(1),'color',color,'vis',vis{1},...
                    'linestyle','none','marker','.','markersize',7,'erasemode','none');
                if draw2
                    linh(j,2) = line(NaN,NaN,'parent',ax(2),'color',color,'vis',vis{2},...
                        'linestyle','none','marker','.','markersize',7,'erasemode','none');
                end
            else
                linh(j,1) = animatedline('parent',ax(1),'color',color,'vis',vis{1},...
                    'linestyle','none','marker','.','markersize',7);
                if draw2
                    linh(j,2) = animatedline('parent',ax(2),'color',color,'vis',vis{2},...
                        'linestyle','none','marker','.','markersize',7);
                end
            end
            
            % Adaptively refine to make smooth curve
            iter = 0;
            while (any(new)) && (iter < maxrefn)
                drawnow
                neww = extermap.demap(zp(new),w,beta,z,c,qdat);
                wp(new) = neww;
                iter = iter + 1;
                
                % Update the points to show progress
                if verLessThan('matlab','8.4')
                    set(linh(j,1),'xdata',real(wp),'ydata',imag(wp))
                    if draw2
                        set(linh(j,2),'xdata',real(zp),'ydata',imag(zp))
                    end
                else
                    addpoints(linh(j,1),real(wp(new)),imag(wp(new)))
                    if draw2
                        addpoints(linh(j,1),real(zp(new)),imag(zp(new)))
                    end
                end
                drawnow
                
                % Add points to zp where necessary
                [zp,wp,new] = scpadapt(zp,wp,minlen,maxlen,axis);
            end
            
            % Set the lines to be solid
            if verLessThan('matlab','8.4')
                set(linh(j,1),'erasemode','back')
                set(linh(j,1),'marker','none','linestyle','-','user',zp)
                if draw2
                    % Replace the points with just the ends
                    set(linh(j,2),'erasemode','back')
                    set(linh(j,2),'marker','none','linestyle','-',...
                        'xdata',[0 1]*cos(theta(j)),'ydata',[0 1]*sin(theta(j)))
                end
            else
                clearpoints(linh(j,1))
                addpoints(linh(j,1),real(wp),imag(wp));
                set(linh(j,1),'marker','none','linestyle','-','user',zp)
                if draw2
                    % Replace the points with just the ends
                    clearpoints(linh(j,2))
                    addpoints(linh(j,2),[0 1]*cos(theta(j)),[0 1]*sin(theta(j)))
                    set(linh(j,2),'marker','none','linestyle','-')
                end
            end
            drawnow
        end
        
        linh = [linh1;linh];
        if ~draw2
            linh = linh(:,1);
        end
        if verLessThan('matlab','8.4'), set(linh,'erasemode','normal'),end
        
        drawnow
        
        if turn_off_hold, hold off, end;
        if nargout > 0
            H = linh;
            if nargout > 1
                R2 = R;
                if nargout > 2
                    THETA = theta;
                end
            end
        end
    end
    
    function I = dequad(z1,z2,sing1,z,beta,qdat)
        %DEQUAD (not intended for calling directly by the user)
        %   Numerical quadrature for the exterior map.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: dequad.m 212 2002-09-25 17:31:37Z driscoll $
        
        %   z1,z2 are vectors of left and right endpoints.  sing1 is a vector of
        %   integer indices which label the singularities in z1.  So if sing1(5)
        %   = 3, then z1(5) = z(3).  A zero means no singularity.  z is the
        %   vector of prevertices (all singularities except the origin); beta is
        %   the vector of associated turning angles.  qdat is quadrature data
        %   from SCQDATA.
        %
        %   Make sure that z and beta are column vectors.
        %
        %   DEQUAD integrates from a possible singularity at the left end to a
        %   regular point at the right.  If both endpoints are singularities,
        %   you must break the integral into two pieces and make two calls.
        %
        %   The integral is subdivided, if necessary, so that no singularity
        %   lies closer to the left endpoint than 1/2 the length of the
        %   integration (sub)interval.  But the singularity at the origin is NOT
        %   accounted for in this decision.
        
        nqpts = size(qdat,1);
        n = length(z);
        bigz = z(:,ones(1,nqpts));
        beta = [beta(:);-2];
        bigbeta = beta(:,ones(1,nqpts));
        if isempty(sing1)
            sing1 = zeros(length(z1),1);
        end
        
        I = zeros(size(z1));
        nontriv = find(z1(:)~=z2(:))';
        
        for k = nontriv
            za = z1(k);
            zb = z2(k);
            sng = sing1(k);
            
            % Allowable integration step, based on nearest singularity.
            dist = min(1,2*min(abs(z([1:sng-1,sng+1:n])-za))/abs(zb-za));
            zr = za + dist*(zb-za);
            % Adjust Gauss-Jacobi nodes and weights to interval.
            ind = sng + (n+1)*(sng==0);
            nd = ((zr-za)*qdat(:,ind) + zr + za).'/2; % nodes
            wt = ((zr-za)/2) * qdat(:,ind+n+1);	% weights
            terms = 1 - nd(ones(n,1),:)./bigz;
            if any(terms(:)==0)
                % Endpoints are practically coincident.
                I(k) = 0;
            else
                terms = [terms;nd];
                % Use Gauss-Jacobi on first subinterval, if necessary.
                if sng > 0
                    terms(sng,:) = terms(sng,:)./abs(terms(sng,:));
                    wt = wt*(abs(zr-za)/2)^beta(sng);
                end
                I(k) = exp(sum(log(terms).*bigbeta))*wt;
                while dist < 1
                    % Do regular Gaussian quad on other subintervals.
                    zl = zr;
                    dist = min(1,2*min(abs(z-zl))/abs(zl-zb));
                    zr = zl + dist*(zb-zl);
                    nd = ((zr-zl)*qdat(:,n+1) + zr + zl).'/2;
                    wt = ((zr-zl)/2) * qdat(:,2*n+2);
                    terms = 1 - nd(ones(n,1),:)./bigz;
                    terms = [terms;nd];
                    I(k) = I(k) + exp(sum(log(terms).*bigbeta))*wt;
                end
            end
        end
    end
end

end
