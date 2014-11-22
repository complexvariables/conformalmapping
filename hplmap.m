classdef hplmap < scmap
%HPLMAP Schwarz-Christoffel half-plane map object.
%   HPLMAP(P) constructs a Schwarz-Christoffel half-plane map object for
%   the polygon P. The parameter problem is solved using default options
%   for the prevertices and the multiplicative constant.
%
%   HPLMAP(P,OPTIONS) uses an options structure of the type created by
%   SCMAPOPT in solving the parameter problem.
%
%   HPLMAP(P,Z) creates a hplmap object having the given prevertices Z
%   (the mulitiplicative constant is found automatically).
%   HPLMAP(P,Z,C) also uses the given constant. An OPTIONS argument can
%   be added, although only the error tolerance will be used.
%
%   HPLMAP(M), where M is a hplmap object, just returns M.
%
%   HPLMAP(M,P) returns a new hplmap object for the polygon P using the
%   options in hplmap M. The prevertices of M will be used as the
%   starting guess for the parameter problem of the new map. Thus P
%   should properly be a perturbation of the polygon for M. An OPTIONS
%   structure may also be given, to override the options in M.
%
%   HPLMAP(Z,ALPHA) creates a map using the given prevertices and the
%   interior polygon angles described by ALPHA (see POLYGON help). The
%   image polygon is deduced by computing S-C integrals assuming a
%   multiplicative constant of 1. HPLMAP(Z,ALPHA,C) uses the given
%   constant instead.
%
%   See also SCMAPOPT, classes POLYGON, SCMAP.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

%   Copyright 1998-2001 by Toby Driscoll.
%   $Id: hplmap.m 215 2002-10-23 18:19:50Z driscoll $    

properties
    prevertex
    constant
    qdata
    accuracyVar
end

methods
    function map = hplmap(varargin)
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
            case 'hplmap'
                oldmap = varargin{1};
                % Continuation of given map to given polygon
                poly = varargin{2};
                opt = scmapopt(oldmap);
                z0 = prevertex(oldmap);
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
                        z = z(:);
                    elseif length(arg) == 1
                        c = arg;
                    else
                        msg = 'Unable to parse argument ''%s''.';
                        error(sprintf(msg,inputname(j+1)))
                    end
                end
                
            case 'double'
                % Args are the prevertex vector, then angle vector
                z = varargin{1}(:);
                alpha = varargin{2}(:);
                if ~isinf(z(end))
                    z = [z;Inf];
                    alpha = [alpha;1];
                end
                poly = polygon(NaN*alpha*1i,alpha);  %  nonsense vertices
                c = 1;
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
                msg = 'Expected ''%s'' to be a polygon, hplmap, or prevertex vector.';
                error(msg,inputname(1))
                
        end % switch
        
        
        % Retrieve options
        opt = scmapopt(opt);
        
        % Take actions based on what needs to be filled in
        
        if isempty(z)
            [w,beta] = scfix('hp',vertex(poly),angle(poly)-1);
            poly = polygon(w,beta+1);
            
            [z,c,qdata] = hplmap.hpparam(w,beta,opt.InitialGuess,opt);
        end
        
        if isempty(qdata)
            % Base accuracy of quadrature on given options
            nqpts = ceil(-log10(opt.Tolerance));
            alpha = angle(poly);
            qdata = scqdata(alpha(1:end-1)-1,nqpts);
        end
        
        if isempty(c)
            % Find constant
            w = vertex(poly);
            beta = angle(poly)-1;
            idx = 1 + find(~isinf(z(2:end)), 1 );
            mid = mean(z([1 idx])) + 1i*diff(real(z([1 idx])))/2;
            I = hplmap.hpquad(z(1),mid,1,z(1:end-1),beta(1:end-1),qdata) - ...
                hplmap.hpquad(z(idx),mid,idx,z(1:end-1),beta(1:end-1),qdata);
            c = diff(w([1 idx]))/I;
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
        map.accuracyVar = accuracy(map);        
    end
    
    function acc = accuracy(M)
        %ACCURACY Apparent accuracy of Schwarz-Christoffel half-plane map.
        %   ACCURACY(M) estimates the accuracy of the Schwarz-Christoffel
        %   half-plane map M. The technique used is to compare the differences
        %   between successive finite vertices to the integral between the
        %   corresponding prevertices, and return the maximum.
        %
        %   See also HPLMAP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: accuracy.m 7 1998-05-10 04:37:19Z tad $
        
        % If an accuracy has been assigned, don't question it
        if ~isempty(M.accuracyVar)
            acc = M.accuracyVar;
            return
        end
        
        % Get data for low-level functions
        p = M.polygon;
        w = vertex(p);
        n = length(w);
        beta = angle(p) - 1;
        z = M.prevertex;
        c = M.constant;
        qdata = M.qdata;
        
        % Test accuracy by integrating between consecutive finite prevertices, and
        % comparing to differences of vertices.
        n = length(w);
        idx = find(~isinf(w(1:n-1)));		% exclude last prevert, at Inf
        
        % Two columns hold endpoint indices for integrations
        idx = [idx(1:end-1) idx(2:end)];
        
        % Find midpoints. Go into upper half-plane to avoid integrating through
        % skipped prevertices.
        zz = z(idx).';
        mid = mean(zz).';
        mid = mid + i*abs(diff(zz).')/2;
        
        % Do the integrations
        I = M.hpquad(z(idx(:,1)),mid,idx(:,1),z(1:n-1),beta(1:n-1),qdata) - ...
            M.hpquad(z(idx(:,2)),mid,idx(:,2),z(1:n-1),beta(1:n-1),qdata);
        
        acc = max(abs( c*I - diff(w(idx).').' ));
    end
    
    function out = char(f)
        %CHAR   Pretty-print a Schwarz-Christoffel half-plane map.
        
        %   Copyright 2001 by Toby Driscoll.
        %   $Id: char.m 158 2001-07-20 14:05:59Z driscoll $
        
        p = f.polygon;
        w = vertex(p);
        alpha = angle(p);
        z = f.prevertex;
        c = f.constant;
        
        if length(z) < length(w)
            z = [z(:);Inf];
        end
        
        L = cell(2,1);
        L{1} = '      vertex               alpha         prevertex       ';
        L{2} = ' --------------------------------------------------------';
        u = real(w);
        v = imag(w);
        fmt = ' %8.5f %c %7.5fi     %8.5f    %20.12e';
        for j = 1:length(w)
            if v(j) < 0
                s = '-';
            else
                s = '+';
            end
            L{end+1} = sprintf(fmt,u(j),s,abs(v(j)),alpha(j),z(j));
        end
        L{end+1} = ' ';
        if imag(c) < 0
            s = '-';
        else
            s = '+';
        end
        L{end+1} = sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c)));
        L{end+1} = sprintf('  Apparent accuracy is %.2e',f.accuracyVar);
        L{end+1} = ' ';
        
        out = L;
    end
    
    function M1 = diskmap(M)
        %DISKMAP Convert Schwarz-Christoffel half-plane map to a map from the disk.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: diskmap.m 7 1998-05-10 04:37:19Z tad $
        
        p = M.polygon;
        [z1,c1] = M.hp2disk(vertex(p),angle(p)-1,M.prevertex,M.constant);
        M1 = diskmap(p,scmapopt(M),z1,c1);
    end
    
    function out = display(M)
        %DISPLAY Display parameters of a Schwarz-Christoffel half-plane map.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: display.m 118 2001-05-03 21:18:27Z driscoll $
        
        p = M.polygon;
        w = vertex(p);
        alpha = angle(p);
        z = M.prevertex;
        c = M.constant;
        
        if length(z) < length(w)
            z = [z(:);Inf];
        end
        
        L = {' '; '  hplmap object:'; ' '};
        L{4} = '      vertex               alpha         prevertex       ';
        L{5} = ' --------------------------------------------------------';
        u = real(w);
        v = imag(w);
        fmt = ' %8.5f %c %7.5fi     %8.5f    %20.12e';
        for j = 1:length(w)
            if v(j) < 0
                s = '-';
            else
                s = '+';
            end
            L{end+1} = sprintf(fmt,u(j),s,abs(v(j)),alpha(j),z(j));
        end
        L{end+1} = ' ';
        if imag(c) < 0
            s = '-';
        else
            s = '+';
        end
        L{end+1} = sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c)));
        L{end+1} = ' ';
        L{end+1} = sprintf('  Apparent accuracy is %.2e',M.accuracyVar);
        L{end+1} = ' ';
        
        
        if nargout==0
            fprintf('%s\n',L{:})
        else
            out = L;
        end
    end
    
    function wp = eval(M,zp,tol)
        %EVAL Evaluate Schwarz-Christoffel half-plane map at points.
        %   EVAL(M,ZP) evaluates the Schwarz-Christoffel map M at the points
        %   ZP in the upper half-plane. The default tolerance of M is used.
        %
        %   EVAL(M,ZP,TOL) attempts to give an answer accurate to TOL. If TOL
        %   is less than the accuracy of M, this is unlikely to be met.
        %
        %   See also HPLMAP, EVALINV.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: eval.m 7 1998-05-10 04:37:19Z tad $
        
        p = M.polygon;
        n = length(p);
        
        if nargin < 3
            qdata = M.qdata;
        else
            qdata = tol;
        end
        
        wp = NaN*zp;
        idx = imag(zp) > -eps;
        wp(idx) = M.hpmap(zp(idx),vertex(p),angle(p)-1,M.prevertex,M.constant,qdata);
    end
    
    function fp = evaldiff(M,zp)
        %EVALDIFF Derivative of Schwarz-Christoffel half-plane map at points.
        %   EVALDIFF(M,ZP) computes the derivative of the Schwarz-Christoffel
        %   half-plane map M at the points ZP.
        %
        %   See also HPLMAP, EVAL.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: evaldiff.m 7 1998-05-10 04:37:19Z tad $
        
        z = M.prevertex;
        c = M.constant;
        beta = angle(M.polygon) - 1;
        
        fp = M.hpderiv(zp,z,beta,c);
    end
    
    function zp = evalinv(M,wp,tol,z0)
        %EVALINV Invert Schwarz-Christoffel half-plane map at points.
        %   EVALINV(M,WP) evaluates the inverse of the Schwarz-Christoffel
        %   map M at the points WP in the polygon. The default tolerance of M
        %   is used.
        %
        %   EVALINV(M,WP,TOL) attempts to give an answer accurate to TOL. If
        %   TOL is smaller than the accuracy of M, this is unlikely to be met.
        %
        %   EVALINV(M,WP,TOL,Z0) uses given starting points. Z0 must be either
        %   the same size as WP or a complex scalar (to be expanded to that
        %   size). It is used for the starting approximation to the inverse
        %   image of WP. The starting guess need not be close to the correct
        %   answer; however, the straight line segment between WP(K) and the
        %   forward image of Z0(K) must lie entirely inside the polygon, for
        %   each K.
        %
        %   See also HPLMAP, EVAL.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: evalinv.m 91 2000-05-17 23:05:51Z tad $
        
        % Assign empties to missing input args
        import sctool.*
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
            tol = M.accuracyVar;
        end
        
        if ~isempty(z0)
            if length(z0) == 1
                %z0 = repmat(z0,size(wp));
            elseif any(size(z0) ~= size(wp))
                msg = 'Argument %s must be a complex scalar or the same size as %s.';
                error(sprintf(msg,inputname(z0),inputname(1)));
            end
        end
        
        p = M.polygon;
        n = length(p);
        w = vertex(p);
        beta = angle(p) - 1;
        z = M.prevertex;
        c = M.constant;
        
        zp = NaN*wp;
        %idx = logical(isinpoly(wp,p));
        idx = logical(ones(size(wp)));
        zp(idx) = M.hpinvmap(wp(idx),w,beta,z,c,qdata,z0,[0 tol]);
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
        %   Given an hplmap M, FORWARDPOLY(M) returns the polygon that is
        %   formed using the prevertices, angles, and quadrature data of that
        %   map. If the prevertices were found from the solution of a
        %   parameter problem, then the result should agree closely with the
        %   original polygon that was supplied.
        
        %   Copyright (c) 1998 by Toby Driscoll.
        %   $Id: forwardpoly.m 37 1998-06-30 00:22:40Z tad $
        
        z = map.prevertex;
        alpha = angle(map.polygon);
        c = map.constant;
        
        n = length(z);
        
        % Since there is no parameter problem, use high accuracy in quadrature.
        qdata = scqdata(alpha(1:n-1)-1,16);
        
        w = zeros(n,1);
        atinf = (alpha < eps);
        w(atinf) = Inf;
        
        % Endpoints of integrations. Because the last prevertex is at Inf, we
        % shouldn't try to integrate there.
        idx = find(~atinf);
        if idx(end)==n, idx(end) = []; end
        endpt = [idx(1:end-1) idx(2:end)];
        
        % Midpoints are in upper half-plane. Always make 45 degrees with real line.
        mid = mean(z(endpt),2) + i*diff(z(endpt),1,2)/2;
        
        % Integrations
        I = map.hpquad(z(endpt(:,1)),mid,endpt(:,1),z(1:n-1),alpha(1:n-1)-1,qdata) - ...
            map.hpquad(z(endpt(:,2)),mid,endpt(:,2),z(1:n-1),alpha(1:n-1)-1,qdata);
        
        % Deduce vertices
        w(idx) = c*cumsum([0;I]);
        
        % Get the last vertex via intersection
        if alpha(n) > 0
            if abs(alpha(n)-1) < 5*eps | abs(alpha(n)-2) < 5*eps
                error(['Cannot deduce last vertex when its adjacent sides are' ...
                    ' collinear.'])
            elseif any(atinf([1 2 n-1]))
                error('Vertices 1, 2, and end-1 must be finite.')
            else
                % Here's the direction from w(1)
                d1 = (w(2)-w(1))*exp(i*pi*alpha(1));
                % Get the direction from w(n-1)
                d2 = angle(w(2)-w(1)) + sum(pi*(1-alpha(2:n-1)));
                d2 = exp(i*d2);
                b = w(n-1) - w(1);
                s = [real([d1 -d2]);imag([d1 -d2])]\[real(b);imag(b)];
                w(n) = w(1) + s(1)*d1;
            end
        end
        
        p = polygon(w,alpha);
    end
    
    function varargout = get(map,varargin)
        %GET    Get map parameters.
        %   [VAL1,VAL2,...] = GET(F,'PROP1','PROP2',...) returns the values of the
        %   map F corresponding to the requested properties. Valid properties
        %   are:
        %
        %       polygon, options, prevertex, constant
        
        % Copyright 1999-2003 by Toby Driscoll.
        % $Id: get.m 237 2003-01-15 15:29:15Z driscoll $
        
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
        if isa(M,'double') & isa(c,'hplmap')
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
        
        v.prevertex = M.prevertex;
        v.constant = M.constant;
    end
    
    function [h,re,im] = plot(M,varargin)
        %PLOT Visualize a Schwarz-Christoffel half-plane map.
        %   PLOT(M) plots the polygon associated with the Schwarz-Christoffel
        %   half-plane map M and the images of ten evenly spaced vertical rays
        %   and horizontal lines under the S-C transformation.
        %
        %   PLOT(M,NRE,NIM) plots the images of NRE vertical rays and NIM
        %   horizontal lines.
        %
        %   PLOT(M,RE,IM) plots the vertical rays at abscissae given by the
        %   entries of RE and horizontal lines at the ordinates specified in IM.
        %
        %   PLOT(M,TOL) or PLOT(M,NRE,NIM,TOL) or PLOT(M,RE,IM,TOL)
        %   computes the map with accuracy roughly TOL. Normally TOL defaults to
        %   1e-4 or the accuracy of M, whichever is greater.
        %
        %   See also HPLMAP, EVAL.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: plot.m 7 1998-05-10 04:37:19Z tad $
        
        p = M.polygon;
        w = vertex(p);
        beta = angle(p) - 1;
        z = M.prevertex;
        c = M.constant;
        
        if nargin == 1
            [a1,a2,a3] = M.hpplot(w,beta,z,c);
        elseif length(varargin) == 1
            % Tolerance given only
            [a1,a2,a3] = M.hpplot(w,beta,z,c,10,10,ceil(-log10(varargin{1})));
        elseif length(varargin) == 2
            % RE,IM given only
            [a1,a2,a3] = M.hpplot(w,beta,z,c,varargin{1},varargin{2});
        else
            % All given
            nqpts = ceil(-log10(varargin{3}));
            [a1,a2,a3] = M.hpplot(w,beta,z,c,varargin{1},varargin{2},nqpts);
        end
        
        if nargout > 0
            h = a1;
            re = a2;
            im = a3;
        end
    end
end

methods(Hidden, Static)
    function [zd,cd] = hp2disk(w,beta,z,c)
        %HP2DISK Convert solution from the half-plane to one from the disk.
        %   [Z,C] = HP2DISK(W,BETA,Z,C) quickly transforms the solution Z,C of
        %   the Schwarz-Christoffel half-plane mapping parameter problem to the
        %   solution ZD,CD of the disk problem.
        %
        %   See also DISK2HP, HPPARAM, DPARAM.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: hp2disk.m 298 2009-09-15 14:36:37Z driscoll $
        
        n = length(w);
        zd = zeros(size(z));
        if isinf(z(n))
            zd(n) = 1;
            zd(1:n-1) = (z(1:n-1)-i)./(z(1:n-1)+i);
        else
            zd = (z-i)./(z+i);
            zd = zd/zd(n);
        end
        zd = sign(zd);
        
        % Recalculate constant from scratch.
        mid = (zd(1)+zd(2))/2;
        qdat = sctool.scqdata(beta,16);
        cd = (w(1) - w(2))/...
            (diskmap.dquad(zd(2),mid,2,zd,beta,qdat) - diskmap.dquad(zd(1),mid,1,zd,beta,qdat));
    end
    
    function fprime = hpderiv(zp,z,beta,c)
        %HPDERIV Derivative of the half-plane map.
        %   HPDERIV(ZP,Z,BETA,C) returns the derivative at the points of ZP of
        %   the Schwarz-Christoffel half-plane map defined by Z, BETA, and C.
        %
        %   See also HPPARAM, HPMAP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: hpderiv.m 298 2009-09-15 14:36:37Z driscoll $
        
        % Support old syntax
        if nargin < 4
            c = 1;
        end
        
        zf = z(~isinf(z));
        beta = beta(~isinf(z));
        zprow = zp(:).';
        fprime = zeros(size(zp));
        
        npts = length(zp(:));
        terms = zprow(ones(length(beta),1),:) - zf(:,ones(npts,1));
        fprime(:) = c*exp(sum(log(terms).*beta(:,ones(npts,1))));
    end
    
    function hpdisp(w,beta,z,c)
        %HPDISP Display results of Schwarz-Christoffel half-plane parameter problem.
        %   HPDISP(W,BETA,Z,C) displays the results of HPPARAM in a pleasant
        %   way.
        %
        %   See also HPPARAM, HPPLOT.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: hpdisp.m 298 2009-09-15 14:36:37Z driscoll $
        
        if length(z) < length(w)
            z = [z(:);Inf];
        end
        disp(' ')
        disp('      vertex [w]           beta          prevertex [z]   ')
        disp(' --------------------------------------------------------')
        u = real(w);
        v = imag(w);
        for j = 1:length(w)
            if v(j) < 0
                s = '-';
            else
                s = '+';
            end
            disp(sprintf(' %8.5f %c %7.5fi     %8.5f    %20.12e',...
                u(j),s,abs(v(j)),beta(j),z(j)));
        end
        disp(' ')
        if imag(c) < 0
            s = '-';
        else
            s = '+';
        end
        disp(sprintf('  c = %.8g %c %.8gi\n',real(c),s,abs(imag(c))))
    end
    
    function zdot = hpimapfun(wp,yp,flag,scale,z,beta,c)
        %   Used by HPINVMAP for solution of an ODE.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: hpimapfun.m 298 2009-09-15 14:36:37Z driscoll $
        
        lenyp = length(yp);
        lenzp = lenyp/2;
        
        % Don't allow points in lower half-plane. This really messes up the
        % derivative calculation.
        zp = yp(1:lenzp) + i*max(0,yp(lenzp+1:lenyp));
        
        f = scale./hplmap.hpderiv(zp,z,beta,c);
        zdot = [real(f);imag(f)];
    end
    
    function zp = hpinvmap(wp,w,beta,z,c,qdat,z0,options)
        %HPINVMAP Schwarz-Christoffel half-plane inverse map.
        %   HPINVMAP(WP,W,BETA,Z,C,TOL) computes the inverse of the
        %   Schwarz-Christoffel half-plane map (i.e., from the polygon to the
        %   upper half-plane ) at the points given in vector WP. The other
        %   arguments are as in HPPARAM. TOL is a scalar tolerance, or a
        %   quadrature-data matrix QDAT as returned by SCQDATA, or may be
        %   omitted.
        %
        %   The default algorithm is to solve an ODE in order to obtain a fair
        %   approximation for ZP, and then improve ZP with Newton
        %   iterations. The ODE solution at WP requires a vector Z0 whose
        %   forward image W0 is such that for each j, the line segment
        %   connecting WP(j) and W0(j) lies inside the polygon. By default Z0 is
        %   chosen by a fairly robust automatic process. Using a parameter (see
        %   below), you can choose to use either an ODE solution or Newton
        %   iterations exclusively.
        %
        %   HPINVMAP(WP,W,BETA,Z,C,TOL,Z0) has two interpretations. If the ODE
        %   solution is being used, Z0 overrides the automatic selection of
        %   initial points. (This can be handy in convex polygons, where the
        %   choice of Z0 is trivial.) Otherwise, Z0 is taken as an initial guess
        %   to ZP. In either case, if length(Z0)==1, the value Z0 is used for
        %   all elements of WP; otherwise, length(Z0) should equal length(WP).
        %
        %   HPINVMAP(WP,W,BETA,Z,C,TOL,Z0,OPTIONS) uses a vector of parameters
        %   that control the algorithm. See SCINVOPT.
        %
        %   See also SCINVOPT, HPPARAM, HPMAP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: hpinvmap.m 298 2009-09-15 14:36:37Z driscoll $
        
        import sctool.*
        n = length(w);
        w = w(:);
        beta = beta(:);
        z = z(:);
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
        
        z0 = z0(:);
        
        [ode,newton,tol,maxiter] = scinvopt(options);
        
        nfin = n - isinf(z(n));
        if isempty(qdat)
            qdat = tol;
        end
        
        if length(qdat)==1
            qdat = scqdata(beta(1:n-1),max(ceil(-log10(qdat)),2));
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
                mapfun = @(zp) hplmap.hpmap(zp,w,beta,z,c,qdat);
                [z0,w0] = sctool.findz0('hp',wp(~done),mapfun,w,beta,z,c,qdat);
            else
                w0 = hplmap.hpmap(z0,w,beta,z,c,qdat);
                if length(z0)==1 && lenwp > 1
                    z0 = z0(:,ones(lenwp,1)).';
                    w0 = w0(:,ones(lenwp,1)).';
                end
                %%w0 = w0(~done);
                %%z0 = z0(~done);
            end
            
            % Use relaxed ODE tol if improving with Newton.
            odetol = max(tol,1e-3*(newton));
            
            % Rescale dependent coordinate
            scale = (wp(~done) - w0(:));
            
            % Solve ODE
            z0 = [real(z0);imag(z0)];
            odefun = @(w,y) hplmap.hpimapfun(w,y,scale,z,beta,c);
            [t,y] = ode23(odefun,[0,0.5,1],z0,odeset('abstol',odetol));
            [m,leny] = size(y);
            zp(~done) = y(m,1:lenwp)+sqrt(-1)*y(m,lenwp+1:leny);
            out = imag(zp) < 0;
            zp(out) = real(zp(out));
            
            clear global SCIMDATA
        end
        
        % Newton iterations
        if newton
            if ~ode
                zn = z0(:);
                if length(z0)==1 & lenwp > 1
                    zn = zn(:,ones(lenwp,1));
                end
                zn(done) = zp(done);
            else
                zn = zp(:);
            end
            
            wp = wp(:);
            k = 0;
            while ~all(done) & k < maxiter
                F = wp(~done) - hplmap.hpmap(zn(~done),w,beta,z,c,qdat);
                m = length(F);
                dF = hplmap.hpderiv(zn(~done),z,beta,c);
                znew = zn(~done) + F(:)./dF(:);
                zn(~done) = real(znew) + i*max(0,imag(znew));
                done(~done) =  (abs(F) < tol);
                k = k + 1;
            end
            if any(abs(F)> tol)
                str = sprintf('Check solution; maximum residual = %.3g\n',max(abs(F)));
                warning(str)
            end
            zp(:) = zn;
        end
    end

% Previously in hpinvmap.m -- seems to duplicate above of same name. -- EK
%
%     function zdot = hpimapfun(wp,yp,scale,z,beta,c);
%         %   Used by HPINVMAP for solution of an ODE.
%         
%         %   Copyright 1998 by Toby Driscoll.
%         %   $Id: hpimapfun.m 298 2009-09-15 14:36:37Z driscoll $
%         
%         lenyp = length(yp);
%         lenzp = lenyp/2;
%         
%         % Don't allow points in lower half-plane. This really messes up the
%         % derivative calculation.
%         zp = yp(1:lenzp) + i*max(0,yp(lenzp+1:lenyp));
%         
%         f = scale./hpderiv(zp,z,beta,c);
%         zdot = [real(f);imag(f)];
%     end
    
    function wp = hpmap(zp,w,beta,z,c,qdat)
        %HPMAP  Schwarz-Christoffel half-plane map.
        %   HPMAP(ZP,W,BETA,Z,C,QDAT) computes the values of the
        %   Schwarz-Christoffel half-plane map at the points in vector ZP.  The
        %   polygon's vertices should be given in W and the arguments Z, C, and
        %   QDAT should be computed by HPPARAM.  HPMAP returns a vector the same
        %   size as ZP.
        %
        %   HPMAP(ZP,W,BETA,Z,C,TOL) uses quadrature data intended to give an
        %   answer accurate to within TOL.
        %
        %   HPMAP(ZP,W,BETA,Z,C) uses a tolerance of 1e-8.
        %
        %   See also HPPARAM, HPPLOT, HPINVMAP.

        %   Copyright 1998 by Toby Driscoll.
        %   $Id: hpmap.m 298 2009-09-15 14:36:37Z driscoll $

        if isempty(zp)
            wp = [];
            return
        end
        import sctool.*

        n = length(w);
        w = w(:);
        beta = beta(:);
        z = z(:);

        % Quadrature data and error tolerance
        if nargin < 6
            tol = 1e-8;
            qdat = scqdata(beta(1:n-1),8);
        elseif length(qdat)==1
            tol = qdat;
            qdat = scqdata(beta(1:n-1),max(ceil(-log10(tol)),8));
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
        zpinf = isinf(zp);
        wp(zpinf) = w(n);
        vertex = vertex | zpinf;

        % "Bad" points are closest to a prevertex of infinity.
        atinf = find(isinf(w)); 		% infinite vertices
        bad = ismember(sing,atinf) & ~vertex;

        if any(bad)
            % Can't integrate starting at pre-infinity: which neighboring prevertex
            % to use?
            direcn = real(zp(bad)-z(sing(bad)));
            sing(bad) = sing(bad) + sign(direcn) + (direcn==0);
            % Midpoints of these integrations
            mid = (z(sing(bad)) + zp(bad)) / 2;
        end

        % zs = the starting singularities
        zs = z(sing);
        % ws = f(zs)
        ws = w(sing);

        % Compute the map directly at "normal" points.
        normal = ~bad & ~vertex;
        if any(normal)
            I = hplmap.hpquad(zs(normal),zp(normal),sing(normal),z(1:n-1),beta(1:n-1),qdat);
            wp(normal) = ws(normal) + c*I;
        end

        % Compute map at "bad" points, in stages. Stop at midpoint to avoid
        % integration where right endpoint is close to a singularity.
        if any(bad)
            I1 = hplmap.hpquad(zs(bad),mid,sing(bad),z(1:n-1),beta(1:n-1),qdat);
            I2 = -hplmap.hpquad(zp(bad),mid,zeros(sum(bad),1),z(1:n-1),beta(1:n-1),qdat);
            wp(bad) = ws(bad) + c*(I1 + I2);
        end

        wp = reshape(wp,shape);
    end
    
    function [z,c,qdat] = hpparam(w,beta,z0,options)
        %HPPARAM Schwarz-Christoffel half-plane parameter problem.
        %   [Z,C,QDAT] = HPPARAM(W,BETA) solves the Schwarz-Christoffel
        %   parameter problem with the upper half-plane as fundamental domain
        %   and interior of the specified polygon as the target. W must be a
        %   vector of the vertices of the polygon, specified in counterclockwise
        %   order. BETA is a vector of turning angles; see SCANGLES. If
        %   successful, HPPARAM will return Z, a vector of the pre-images of W;
        %   C, the multiplicative constant of the conformal map; and QDAT, an
        %   optional matrix of quadrature data used by some of the other S-C
        %   routines.
        %
        %   [Z,C,QDAT] = HPPARAM(W,BETA,Z0) uses Z0 as an initial guess for Z.
        %
        %   [Z,C,QDAT] = HPPARAM(W,BETA,TOL) attempts to find an answer within
        %   tolerance TOL. (Also see SCPAROPT.)
        %
        %   [Z,C,QDAT] = HPPARAM(W,BETA,Z0,OPTIONS) uses a vector of control
        %   parameters. See SCPAROPT.
        %
        %   See also SCPAROPT, DRAWPOLY, HPDISP, HPPLOT, HPMAP, HPINVMAP.
        
        %   Copyright 1998--2001 by Toby Driscoll.
        %   $Id: hpparam.m 298 2009-09-15 14:36:37Z driscoll $
        
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
        
        err = sccheck('hp',w,beta);
        if err==1
            fprintf('Use SCFIX to make polygon obey requirements\n')
            error(' ')
        end
        
        [trace,tol,method] = parseopt(options);
        if length(z0)==1
            tol = z0;
            z0 = [];
        end
        nqpts = max(ceil(-log10(tol)),4);
        qdat = scqdata(beta(1:n-1),nqpts); 	% quadrature data
        
        atinf = (beta <= -1);
        
        % Find prevertices (solve param problem)
        if n==3
            z = [-1;1;Inf];
            
        else
            
            % Set up normalized lengths for nonlinear equations:
            
            % indices of left and right integration endpoints
            left = find(~atinf(1:n-2));
            right = 1+find(~atinf(2:n-1));
            cmplx = ((right-left) == 2);
            % normalize lengths by w(2)-w(1)
            nmlen = (w(right)-w(left))/(w(2)-w(1));
            % abs value for finite ones; Re/Im for infinite ones
            nmlen = [abs(nmlen(~cmplx));real(nmlen(cmplx));imag(nmlen(cmplx))];
            % first entry is useless (=1)
            nmlen(1) = [];
            
            % Set up initial guess
            if isempty(z0)
                z0 = linspace(-1,1,n-1)';
            else
                z0 = z0(:);
                z0 = z0*2/(z0(n-1)-z0(1));
                z0 = z0-z0(1)-1;
            end
            %y0 = log(diff(z0(1:n-2)));
            y0 = log(diff(z0(1:n-2))./diff(z0(2:n-1)));
            
            % Solve nonlinear system of equations:
            
            % package data
            fdat = {n,beta(1:n-1),nmlen,left,right,logical(cmplx),qdat};
            % set options
            opt = zeros(16,1);
            opt(1) = trace;
            opt(2) = method;
            opt(6) = 100*(n-3);
            opt(8) = tol;
            opt(9) = min(eps^(2/3),tol/10);
            opt(12) = nqpts;
            try
                [y,termcode] = nesolve(@hplmap.hppfun,y0,opt,fdat);
            catch
                % Have to delete the "waitbar" figure if interrupted
                close(findobj(allchild(0),'flat','Tag','TMWWaitbar'));
                error(lasterr)
            end
            if termcode~=1
                warning('Nonlinear equations solver did not terminate normally.')
            end
            
            % Convert y values to z
            %z = [cumsum([-1;exp(y)]);1;Inf];
            cp = cumprod([1;exp(-y)]);
            z = [0;cumsum(cp)] - [flipud(cumsum(flipud(cp)));0];
            z = [z/z(n-1);Inf];
        end
        
        % Determine multiplicative constant
        mid = mean(z(1:2));
        g = hplmap.hpquad(z(2),mid,2,z(1:n-1),beta(1:n-1),qdat) -...
            hplmap.hpquad(z(1),mid,1,z(1:n-1),beta(1:n-1),qdat);
        c = (w(1)-w(2))/g;
        
    end
    
    function F = hppfun(y,fdat)
        %   Returns residual for solution of nonlinear equations.
        %   $Id: hppfun.m 60 1999-01-29 00:49:09Z tad $
        
        [n,beta,nmlen,left,right,cmplx,qdat] = deal(fdat{:});
        
        % Transform y (unconstr. vars) to z (prevertices)
        cp = cumprod([1;exp(-y)]);
        z = [0;cumsum(cp)] - [flipud(cumsum(flipud(cp)));0];
        z = z/z(n-1);
        
        % Compute the integrals
        zleft = z(left);
        zright = z(right);
        mid = mean([zleft.' ; zright.']).';
        % For integrals between non-adjacent singularities, choose intermediate
        % points in the upper half-plane.
        mid(cmplx) = mid(cmplx) + 1i*(zright(cmplx)-zleft(cmplx))/2;
        ints = hplmap.hpquad(zleft,mid,left,z,beta,qdat) - ...
            hplmap.hpquad(zright,mid,right,z,beta,qdat);
        
        if any(ints==0)
            % Singularities were too crowded in practice.
            warning('Severe crowding')
        end
        
        % Compute nonlinear equation residual values.
        F1 = abs(ints(~cmplx));		% F1(1) = abs(ints(1))
        F1 = F1(2:end)/F1(1);
        F2 = ints(cmplx)/ints(1);
        F = [F1;real(F2);imag(F2)] - nmlen;
        
        %%ints = ints/ints(1);
        %%w = [0;cumsum([1;ints(2:end)])];
        %%set(findobj(0,'tag','polydata'),'xd',real(w),'yd',imag(w))
        %%drawnow
        
    end
    
    function [H,RE,IM] = hpplot(w,beta,z,c,re,im,options)
        %HPPLOT Image of cartesian grid under Schwarz-Christoffel half-plane map.
        %   HPPLOT(W,BETA,Z,C) will adaptively plot the images under the
        %   Schwarz-Christoffel exterior map of ten evenly spaced horizontal and
        %   vertical lines in the upper half-plane. The abscissae of the
        %   vertical lines will bracket the finite extremes of Z.  The arguments
        %   are as in HPPARAM.
        %
        %   HPPLOT(W,BETA,Z,C,M,N) will plot images of M evenly spaced vertical
        %   and N evenly spaced horizontal lines.  The spacing will be the same
        %   in both directions.
        %
        %   HPPLOT(W,BETA,Z,C,RE,IM) will plot images of vertical lines whose
        %   real parts are given in RE and horizontal lines whose imaginary
        %   parts are given in IM.  Either argument may be empty.
        %
        %   HPPLOT(W,BETA,Z,C,RE,IM,OPTIONS) allows customization of HPPLOT's
        %   behavior.  See SCPLTOPT.
        %
        %   H = HPPLOT(W,BETA,Z,C,...) returns a vector of handles to all the
        %   curves drawn in the interior of the polygon.  [H,RE,IM] =
        %   HPPLOT(W,BETA,Z,C,...) also returns the abscissae and ordinates of
        %   the lines comprising the grid.
        %
        %   See also SCPLTOPT, HPPARAM, HPMAP, HPDISP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: hpplot.m 298 2009-09-15 14:36:37Z driscoll $
        
        n = length(w);
        w = w(:);
        beta = beta(:);
        z = z(:);
        import sctool.*
        
        % Parse input
        if nargin < 7
            options = [];
            if nargin < 6
                im = [];
                if nargin < 5
                    re = [];
                end
            end
        end
        
        % Empty arguments default to 10
        if isempty([re(:);im(:)])
            re = 10;
            im = 10;
        end
        
        % Integer arguments must be converted to specific values
        if (length(re)==1) && (re == round(re))
            if re < 1
                re = [];
            elseif re < 2
                re = mean(z([1,n-1]));
            else
                m = re;
                re = linspace(z(1),z(n-1),m);
                dre = diff(re(1:2));
                re = linspace(z(1)-dre,z(n-1)+dre,m);
            end
        end
        if (length(im)==1) && (im == round(im))
            if length(re) < 2
                im = linspace(0,4,im+1);
                im(1) = [];
            else
                im = mean(diff(re))*(1:im);
            end
        end
        
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
        
        % Drawing parameters
        [nqpts,minlen,maxlen,maxrefn] = scpltopt(options);
        len = max(diff(get(ax(1),'xlim')),diff(get(ax(1),'ylim')));
        minlen = len*minlen;
        maxlen = len*maxlen;
        if any(isinf(z))
            qdat = scqdata(beta(1:n-1),nqpts);
        else
            qdat = scqdata(beta,nqpts);
        end
        axlim = axis;
        
        color = 'k';
        
        % Plot vertical lines...
        y2 = max(z(n-1),10);
        linh = gobjects(length(re),2);
        for j = 1:length(re)
            % Start evenly spaced
            zp = re(j) + 1i*[linspace(0,y2,15) Inf].';
            new = true(size(zp));
            new(end) = false;
            wp = NaN(length(zp),1);
            wp(end) = w(n);
            
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
            
            % Adaptive refinement to make smooth curve
            iter = 0;
            while (any(new)) && (iter < maxrefn)
                drawnow
                neww = hplmap.hpmap(zp(new),w,beta,z,c,qdat);
                wp(new) = neww;
                iter = iter + 1;
                
                % Update the points to show progress
                if verLessThan('matlab','8.4')
                    set(linh(j,1),'xdata',real(wp),'ydata',imag(wp))
                    if draw2
                        set(linh(j,2),'xdata',real(zp),'ydata',imag(zp))
                    end
                else
                    % On the first pass, include the image of Inf
                    if iter==1
                        addpoints(linh(j,1),real(wp),imag(wp))
                    else
                        addpoints(linh(j,1),real(wp(new)),imag(wp(new)))
                    end
                    if draw2
                        addpoints(linh(j,1),real(zp(new)),imag(zp(new)))
                    end
                end
                drawnow
                
                % Add points to zp where necessary
                [zp,wp,new] = scpadapt(zp,wp,minlen,maxlen,axlim);
            end
            
            % Set the lines to be solid
            if verLessThan('matlab','8.4')
                set(linh(j,1),'erasemode','back')
                set(linh(j,1),'marker','none','linestyle','-','user',zp)
                if draw2
                    % Replace the points with the endpoints
                    set(linh(j,2),'erasemode','back')
                    set(linh(j,2),'marker','none','linestyle','-',...
                        'xdata',re(j)*[1 1],'ydata',[0 imag(zp(end-1))])
                end
            else
                clearpoints(linh(j,1))
                addpoints(linh(j,1),real(wp(~isnan(wp))),imag(wp(~isnan(wp))));
                set(linh(j,1),'marker','none','linestyle','-','user',zp)
                if draw2
                    % Replace the points with (hopefully) a smooth circle
                    clearpoints(linh(j,2))
                    addpoints(linh(j,2),re(j)*[1 1],[0 imag(zp(end-1))])
                    set(linh(j,2),'marker','none','linestyle','-')
                end
                
            end
            drawnow
        end
        
        z1 = min(-10,z(1));
        z2 = max(40,z(n-1));
        linh1 = linh;
        linh = gobjects(length(im),2);
        for j = 1:length(im)
            % Start evenly spaced
            zp = [-Inf linspace(z1,z2,15) Inf].' + 1i*im(j);
            new = true(size(zp));
            new([1 end]) = 0;
            wp = NaN(length(zp),1);
            wp([1 end]) = w(n);
            
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
            
            % Adaptive refinement to make smooth curve
            iter = 0;
            while (any(new)) && (iter < maxrefn)
                drawnow
                neww = hplmap.hpmap(zp(new),w,beta,z,c,qdat);
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
                [zp,wp,new] = scpadapt(zp,wp,minlen,maxlen,axlim);
            end
            
            % Set the lines to be solid
            if verLessThan('matlab','8.4')
                set(linh(j,1),'erasemode','back')
                set(linh(j,1),'marker','none','linestyle','-','user',zp)
                if draw2
                    % Replace the points with just the ends
                    set(linh(j,2),'erasemode','back')
                    set(linh(j,2),'marker','none','linestyle','-',...
                        'xdata',real(zp([2 end-1])),'ydata',im(j)*[1 1])
                end
            else
                clearpoints(linh(j,1))
                addpoints(linh(j,1),real(wp),imag(wp));
                set(linh(j,1),'marker','none','linestyle','-','user',zp)
                if draw2
                    % Replace the points with just the ends
                    clearpoints(linh(j,2))
                    addpoints(linh(j,2),real(zp([2 end-1])),im(j)*[1 1])
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
                RE = re;
                if nargout > 2
                    IM = im;
                end
            end
        end
        
    end
    
    function I = hpquad(z1,z2,varargin)
        %HPQUAD (not intended for calling directly by the user)
        %   Numerical quadrature for the half-plane map.

        %   Copyright 1998 by Toby Driscoll.
        %   $Id: hpquad.m 298 2009-09-15 14:36:37Z driscoll $

        %   HPQUAD(z1,z2,sing1,z,beta,qdat)
        %   z1,z2 are vectors of left and right endpoints.  sing1 is a vector of
        %   integer indices which label the singularities in z1.  So if sing1(5)
        %   = 3, then z1(5) = z(3).  A zero means no singularity.  z is the
        %   vector of finite singularities; beta is the vector of associated
        %   turning angles.  qdat is quadrature data from SCQDATA.
        %
        %   Make sure z and beta are column vectors.
        %
        %   HPQUAD integrates from a possible singularity at the left end to a
        %   regular point at the right.  If both endpoints are singularities,
        %   you must break the integral into two pieces and make two calls, or
        %   call HPQUAD(z1,z2,sing1,sing2,z,beta,qdat) and accept an automatic
        %   choice.
        %
        %   The integral is subdivided, if necessary, so that no singularity
        %   lies closer to the left endpoint than 1/2 the length of the
        %   integration (sub)interval.

        if nargin==7
            % Break into two pieces with recursive call.
            [sing1,sing2,z,beta,qdat] = deal(varargin{:});
            mid = (z1+z2)/2;
            mid = mid + 1i*abs(mid);
            I1 = hplmap.hpquad(z1,mid,sing1,z,beta,qdat);
            I2 = hplmap.hpquad(z2,mid,sing2,z,beta,qdat);
            I = I1-I2;
            return
        else
            [sing1,z,beta,qdat] = deal(varargin{:});
        end

        nqpts = size(qdat,1);
        % Note: Here n is the total number of *finite* singularities; i.e., the
        % number of terms in the product appearing in the integrand.
        n = length(z);
        bigz = z(:,ones(1,nqpts));
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
            % %  if isempty(dist), dist=1; end
            zr = za + dist*(zb-za);
            ind = rem(sng+n,n+1)+1;
            % Adjust Gauss-Jacobi nodes and weights to interval.
            nd = ((zr-za)*qdat(:,ind) + zr + za).'/2; % G-J nodes
            wt = ((zr-za)/2) * qdat(:,ind+n+1); 	% G-J weights
            terms = nd(ones(n,1),:) - bigz;
            if any(terms(:)==0)
                % Endpoints are practically coincident.
                I(k) = 0;
            else
                % Use Gauss-Jacobi on first subinterval, if necessary.
                if sng > 0
                    terms(sng,:) = terms(sng,:)./abs(terms(sng,:));
                    wt = wt*(abs(zr-za)/2)^beta(sng);
                end
                I(k) = exp(sum(log(terms).*bigbeta,1))*wt;
                while dist < 1
                    % Do regular Gaussian quad on other subintervals.
                    zl = zr;
                    dist = min(1,2*min(abs(z-zl))/abs(zl-zb));
                    zr = zl + dist*(zb-zl);
                    nd = ((zr-zl)*qdat(:,n+1) + zr + zl).'/2;
                    wt = ((zr-zl)/2) * qdat(:,2*n+2);
                    terms = nd(ones(n,1),:) - bigz;
                    I(k) = I(k) + exp(sum(log(terms).*bigbeta,1)) * wt;
                end
            end
        end
    end
end

end
