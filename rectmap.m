classdef rectmap < scmap
%RECTMAP Schwarz-Christoffel rectangle map object.
%   RECTMAP(P,CORNER) constructs a Schwarz-Christoffel rectangle map
%   object for the polygon P. CORNER is a four-vector containing the
%   indices of the vertices that are the images of the rectangle's
%   corners. These indices should be specified in counterclockwise
%   order, and the first two should be the endpoints of a long side of
%   the rectangle.
%
%   RECTMAP(P) requires you to choose the corners graphically.
%
%   RECTMAP(P,CORNER,OPTIONS) or RECTMAP(P,OPTIONS) uses an options
%   structure of the type created by SCMAPOPT in solving the parameter
%   problem.
%
%   RECTMAP(P,Z,C,L) constructs a rectmap using the explicitly given
%   prevertices Z, constant C, and strip length L.
%
%   RECTMAP(M), where M is a rectmap object, just returns M.
%
%   RECTMAP(M,P) returns a new rectmap object for the polygon P using
%   the options in rectmap M. The prevertices of M will be used as the
%   starting guess for the parameter problem of the new map. Thus P
%   should properly be a perturbation (continuation) of the polygon for
%   M. An OPTIONS structure may be added to override options in M. There
%   is no opportunity to change the corner indices.
%
%   See also SCMAPOPT.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

%   Copyright 1998--2001 by Toby Driscoll.

properties
    prevertex
    constant
    qdata
    stripL
    accuracyVal
end

methods
    function map = rectmap(poly, varargin)
        % Assign empties to optional args
        corner = [];
        z = [];
        c = [];
        L = [];
        opt = [];
        import sctool.*
        
        % Branch based on class of first argument
        switch class(poly)
            
            case 'rectmap'
                oldmap = poly;
                % Continuation of given map to given polygon
                poly = varargin{1};
                % Use strip prevertices as initial guess
                %z0 = M.strip_prevertex;
                zr = oldmap.prevertex;
                z0 = rectmap.r2strip(zr,zr,oldmap.stripL);
                if length(z0) ~= length(poly)
                    msg = 'Polygon %s must have the same length as that in %s.';
                    error(msg,inputname(2),inputname(1))
                end
                opt = scmapopt(oldmap);
                if nargin > 2
                    opt = scmapopt(opt,varargin{2});
                end
                opt = scmapopt(opt,'initial',z0);
                corner = corners(oldmap);
                
            case 'polygon'
                % Parse optional arguments
                for j = 1:length(varargin)
                    arg = varargin{j};
                    % Each arg is the corner vector, an options struct, or z/c/L
                    if isa(arg,'struct')
                        opt = arg;
                    elseif (length(arg) == 4) && all(round(arg) == arg)
                        corner = arg;
                    elseif length(arg) == length(poly)
                        % In this case, immediately pick up the two constants as well
                        z = arg;
                        z = z(:);
                        if j~=1 or length(varargin)~=3
                            error('Incorrectly specified prevertices, c, and L.')
                        else
                            c = varargin{2};
                            L = varargin{3};
                            break
                        end
                    else
                        msg = 'Unable to parse argument ''%s''.';
                        error(msg,inputname(j+1))
                    end
                end
                
            otherwise
                msg = 'Expected ''%s'' to be of class polygon or rectmap.';
                error(msg,inputname(1))
                
        end % switch
        
        % Retrieve options
        opt = scmapopt(opt);
        
        % Get data for the low-level functions
        w = vertex(poly);
        n = length(w);
        beta = angle(poly) - 1;
        
        if isempty(z)
            % Request corners
            if isempty(corner)
                msg{1} = 'Select the images of the corners of the rectangle.';
                msg{2} = 'Go in counterclockwise order and select a long rectangle edge first.';
                corner = scselect(w,beta,4,'Select corners',msg);
            end
            
            % Apply SCFIX to enforce solver rules
            [w,beta,corner] = scfix('r',w,beta,corner);
            poly = polygon(w,beta+1);
            
            % Solve parameter problem (always necessary)
            [z,c,L,qdata] = rectmap.rparam(w,beta,corner,opt.InitialGuess,opt);
            
        else
            % Prevertices, etc. given. Renumber to conform
            [w,beta,z] = rectmap.rcorners(w,beta,z);
            nqpts = ceil(-log10(opt.Tolerance));
            qdata = scqdata(beta,nqpts);
            poly = polygon(w,beta+1);
        end
        
        map = map@scmap(poly,opt);
        
        map.prevertex = z;
        map.constant = c;
        map.stripL = L;
        map.qdata = qdata;
        
        % Now fill in apparent accuracy
        map.accuracyVal = accuracy(map);
    end
    
    function acc = accuracy(M)
        %ACCURACY Apparent accuracy of Schwarz-Christoffel rectangle map.
        %   ACCURACY(M) estimates the accuracy of the Schwarz-Christoffel
        %   rectangle map M. The technique used is to compare the differences
        %   between successive finite vertices to the integral between the
        %   corresponding prevertices, and return the maximum.
        %
        %   See also RECTMAP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: accuracy.m 7 1998-05-10 04:37:19Z tad $
        
        % If an accuracy has been assigned, don't question it
        if ~isempty(M.accuracyVal)
            acc = M.accuracyVal;
            return
        end
        
        % Get data for low-level functions
        p = M.polygon;
        w = vertex(p);
        n = length(w);
        beta = angle(p) - 1;
        z = M.prevertex;
        c = M.constant;
        L = M.stripL;
        qdata = M.qdata;
        
        % Renumber to put first corner first
        [w,beta,z,corner] = M.rcorners(w,beta,z);
        
        % Map prevertices to strip
        K = max(real(z));
        Kp = max(imag(z));
        zs = M.r2strip(z,z,L);
        zs = real(zs) + i*round(imag(zs));	% put them *exactly* on edges
        
        % Integrate between consecutive finite pairs on bottom and top
        idxbot = find( ~isinf(w(1:corner(3)-1)) );
        idxtop = corner(3)-1 + find( ~isinf(w(corner(3):n)) );
        
        % Two columns hold endpoint indices for integrations
        idx = [idxbot(1:end-1) idxbot(2:end)];
        idx = [idx ; [idxtop(1:end-1) idxtop(2:end)] ];
        
        % Find midpoints. Go into interior of strip to avoid integrating through
        % skipped prevertices.
        zz = zs(idx).';
        mid = mean(zz).';
        mid = real(mid) + i/2;
        
        % As a final check, integrate once across the strip
        [tmp,k] = min(abs( real(zs(idxtop)-zs(1)) ));
        idx = [idx; [1 idxtop(k)]];
        mid(end+1) = mean(zs(idx(end,:)));
        
        % Add in ends of strip
        ends = find(diff(imag(zs([1:n 1]))));
        zq = [zs(1:ends(1));Inf;zs(ends(1)+1:ends(2));-Inf;zs(ends(2)+1:n)];
        bq = [beta(1:ends(1));0;beta(ends(1)+1:ends(2));0;beta(ends(2)+1:n)];
        wq = [w(1:ends(1));NaN;w(ends(1)+1:ends(2));NaN;w(ends(2)+1:n)];
        % Extend qdat with useless columns at ends
        j = [1:ends(1) n+1 ends(1)+1:ends(2) n+1 ends(2)+1:n n+1];
        qdata = qdata(:,[j j+n+1]);
        
        % Do the integrations
        zleft = zs(idx(:,1));
        zright = zs(idx(:,2));
        idx = idx + (idx > ends(1)) + (idx > ends(2));
        I = stripmap.quad(zleft,mid,idx(:,1),zq,bq,qdata) - ...
            stripmap.quad(zright,mid,idx(:,2),zq,bq,qdata);
        
        acc = max(abs( c*I - diff(wq(idx).').' ));
    end
    
    function out = char(f)
        %CHAR   Pretty-print a Schwarz-Christoffel rectangle map.
        
        %   Copyright 1998-2001 by Toby Driscoll.
        %   $Id: char.m 162 2001-07-20 14:33:00Z driscoll $
        
        p = f.polygon;
        w = vertex(p);
        alpha = angle(p);
        z = f.prevertex;
        c = f.constant;
        
        n = length(w);
        % Deduce corner locations
        left = abs(real(z)-min(real(z))) < eps;
        right = abs(real(z)-max(real(z))) < eps;
        top = abs(imag(z)-max(imag(z))) < eps;
        bot = abs(imag(z)-min(imag(z))) < eps;
        corners = find(left+right+top+bot - 1);
        c1 = find(abs(z-max(real(z))) < eps);
        offset = find(corners==c1);
        corners = corners([offset:4,1:offset-1]);
        rect = z(corners);
        
        L = cell(2,1);
        L{1}=' cnr      vertex              alpha               prevertex       ';
        L{2}=' ------------------------------------------------------------------------';
        
        u = real(w);
        v = imag(w);
        for j = 1:length(w)
            if v(j) < 0
                s = '-';
            else
                s = '+';
            end
            cnr = find(j==corners);
            if isempty(cnr)
                cstr = '    ';
            else
                cstr = sprintf('  %i ',cnr);
            end
            if ~imag(z(j))
                L{end+1}=sprintf('%s %8.5f %c %7.5fi    %8.5f   %16.8e',...
                    cstr,u(j),s,abs(v(j)),alpha(j),z(j));
            else
                L{end+1}=sprintf('%s %8.5f %c %7.5fi    %8.5f   %16.8e + %14.8ei',...
                    cstr,u(j),s,abs(v(j)),alpha(j),real(z(j)),imag(z(j)));
            end
        end
        
        L{end+1} = ' ';
        if imag(c) < 0
            s = '-';
        else
            s = '+';
        end
        L{end+1} = sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c)));
        L{end+1} = sprintf('  Conformal modulus = %.8g',imag(rect(2))/rect(1)/2);
        L{end+1} = sprintf('  Apparent accuracy is %.2e',f.accuracyVal);
        L{end+1} = ' ';
        
        out = L;
    end
    
    function corner = corners(M)
        %CORNERS Indices of rectangle/generalized quadrilateral corners.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: corners.m 7 1998-05-10 04:37:19Z tad $
        
        z = M.prevertex;
        tol = 4*eps;
        
        % Find extent of rectangle
        K = max(real(z));
        Kp = max(imag(z));
        
        % First corner is K + 0i
        dif = repmat(z,1,4) - repmat([K K+i*Kp -K+i*Kp -K],length(z),1);
        [tmp,corner] = min(abs(dif));
        
        corner = corner(:);
    end
    
    function out = display(M)
        %DISPLAY Display parameters of a Schwarz-Christoffel rectangle map.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: display.m 118 2001-05-03 21:18:27Z driscoll $
        
        p = M.polygon;
        w = vertex(p);
        alpha = angle(p);
        z = M.prevertex;
        c = M.constant;
        
        n = length(w);
        % Deduce corner locations
        left = abs(real(z)-min(real(z))) < eps;
        right = abs(real(z)-max(real(z))) < eps;
        top = abs(imag(z)-max(imag(z))) < eps;
        bot = abs(imag(z)-min(imag(z))) < eps;
        corners = find(left+right+top+bot - 1);
        c1 = find(abs(z-max(real(z))) < eps);
        offset = find(corners==c1);
        corners = corners([offset:4,1:offset-1]);
        rect = z(corners);
        
        L = {' '; '  rectmap object:'; ' '};
        L{4}=' cnr      vertex              alpha               prevertex       ';
        L{5}=' ------------------------------------------------------------------------';
        
        u = real(w);
        v = imag(w);
        for j = 1:length(w)
            if v(j) < 0
                s = '-';
            else
                s = '+';
            end
            cnr = find(j==corners);
            if isempty(cnr)
                cstr = '    ';
            else
                cstr = sprintf('  %i ',cnr);
            end
            if ~imag(z(j))
                L{end+1}=sprintf('%s %8.5f %c %7.5fi    %8.5f   %16.8e',...
                    cstr,u(j),s,abs(v(j)),alpha(j),z(j));
            else
                L{end+1}=sprintf('%s %8.5f %c %7.5fi    %8.5f   %16.8e + %14.8ei',...
                    cstr,u(j),s,abs(v(j)),alpha(j),real(z(j)),imag(z(j)));
            end
        end
        
        L{end+1} = ' ';
        if imag(c) < 0
            s = '-';
        else
            s = '+';
        end
        L{end+1} = sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c)));
        L{end+1} = sprintf('  Conformal modulus = %.8g',imag(rect(2))/rect(1)/2);
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
        %EVAL Evaluate Schwarz-Christoffel rectangle map at points.
        %   EVAL(M,ZP) evaluates the Schwarz-Christoffel map M at the points ZP
        %   in the source rectangle of M. The default tolerance of M is used.
        %
        %   EVAL(M,ZP,TOL) attempts to give an answer accurate to TOL. If TOL is
        %   less than the accuracy of M, this is unlikely to be met.
        %
        %   See also RECTMAP, EVALINV.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: eval.m 267 2003-05-08 18:07:51Z driscoll $
        
        p = M.polygon;
        n = length(p);
        z = M.prevertex;
        zr = z(corners(M));
        
        if nargin < 3
            qdata = M.qdata;
        else
            qdata = tol;
        end
        
        wp = NaN*zp;
        opt = scmapopt(M);
        idx = find( isinpoly(zp,polygon(zr),opt.Tolerance) );
        wp(idx) = ...
            M.rmap(zp(idx),vertex(p),angle(p)-1,z,M.constant,M.stripL,qdata);
    end
    
    function fp = evaldiff(M,zp)
        %EVALDIFF Derivative of Schwarz-Christoffel rectangle map at points.
        %   EVALDIFF(M,ZP) computes the derivative of the Schwarz-Christoffel
        %   rectangle map M at the points ZP.
        %
        %   See also RECTMAP, EVAL.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: evaldiff.m 7 1998-05-10 04:37:19Z tad $
        
        z = M.prevertex;
        c = M.constant;
        beta = angle(M.polygon) - 1;
        
        fp = M.rderiv(zp,z,beta,M.constant,M.stripL);
    end
    
    function zp = evalinv(M,wp,tol,z0)
        %EVALINV Invert Schwarz-Christoffel rectangle map at points.
        %   EVALINV(M,WP) evaluates the inverse of the Schwarz-Christoffel
        %   rectangle map M at the points WP in the polygon. The default
        %   tolerance of M is used.
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
        %   See also RECTMAP, EVAL.
        
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
        n = length(p);
        w = vertex(p);
        beta = angle(p) - 1;
        z = M.prevertex;
        c = M.constant;
        L = M.stripL;
        
        zp = NaN*wp;
        %idx = logical(isinpoly(wp,p));
        idx = logical(ones(size(wp)));
        zp(idx) = M.rinvmap(wp(idx),w,beta,z,c,L,qdata,z0,[0 tol]);
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
    
    function varargout = get(map,varargin)
        %GET    Get map parameters.
        %   [VAL1,VAL2,...] = GET(F,'PROP1','PROP2',...) returns the values of the
        %   map F corresponding to the requested properties. Valid properties are:
        %
        %       polygon, options, prevertex, constant, L (stripL)
        
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
                case {'l','str'}
                    varargout{j} = map.stripL;
                otherwise
                    warning(sprintf('Property ''%s'' not recognized.\n',varargin{j}))
                    varargout{j} = [];
            end
        end
    end
    
    function mu = modulus(M)
        %MODULUS Conformal modulus of the generalized quadrilateral.
        %   Returns the conformal modulus of the polygon in the rectmap (the
        %   aspect ratio of the source rectangle).
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: modulus.m 7 1998-05-10 04:37:19Z tad $
        
        z = M.prevertex;
        mu = max(imag(z)) / (2*max(real(z)));
    end
    
    function M = mtimes(M,c)
        %   Scale the image of a map by a complex constant.
        
        %   Copyright (c) 1998 by Toby Driscoll.
        %   $Id: mtimes.m 33 1998-06-29 22:35:40Z tad $
        
        % May need to swap arguments
        if isa(M,'double') & isa(c,'rectmap')
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
        v.stripL = M.stripL;
    end
    
    function [h,re,im] = plot(M,varargin)
        %PLOT Visualize a Schwarz-Christoffel rectangle map.
        %   PLOT(M) plots the polygon associated with the Schwarz-Christoffel
        %   rectangle map M and the images of ten evenly spaced vertical and
        %   horizontal line segments under the S-C transformation.
        %
        %   PLOT(M,NRE,NIM) plots the images of NRE vertical and NIM horizontal
        %   line segments.
        %
        %   PLOT(M,RE,IM) plots the vertical line segments at abscissae given by
        %   the entries of RE and horizontal line segments at the ordinates
        %   specified in IM.
        %
        %   PLOT(M,TOL) or PLOT(M,NRE,NIM,TOL) or PLOT(M,RE,IM,TOL) computes the
        %   map with accuracy roughly TOL. Normally TOL defaults to 1e-4 or the
        %   accuracy of M, whichever is greater.
        %
        %   See also RECTMAP, EVAL.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: plot.m 7 1998-05-10 04:37:19Z tad $
        
        p = M.polygon;
        w = vertex(p);
        beta = angle(p) - 1;
        z = M.prevertex;
        c = M.constant;
        L = M.stripL;
        
        if nargin == 1
            [a1,a2,a3] = M.rplot(w,beta,z,c,L);
        elseif length(varargin) == 1
            % Tolerance given only
            [a1,a2,a3] = M.rplot(w,beta,z,c,L,10,10,ceil(-log10(varargin{1})));
        elseif length(varargin) == 2
            % RE,IM given only
            [a1,a2,a3] = M.rplot(w,beta,z,c,L,varargin{1},varargin{2});
        else
            % All given
            nqpts = ceil(-log10(varargin{3}));
            [a1,a2,a3] = M.rplot(w,beta,z,c,L,varargin{1},varargin{2},nqpts);
        end
        
        if nargout > 0
            h = a1;
            re = a2;
            im = a3;
        end
    end
    
    function zr = rectangle(M)
        %RECTANGLE Return the corners of the rectangle in the fundamental domain.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: rectangle.m 7 1998-05-10 04:37:19Z tad $
        
        zr = prevertex(M);
        zr = zr(corners(M));
    end
end

methods(Hidden,Static)
    function [sn,cn,dn] = ellipjc(u,L,flag)
        %ELLIPJC Jacobi elliptic functions for complex argument.
        %   [SN,CN,DN] = ELLIPJC(U,L) returns the values of the Jacobi
        %   elliptic functions evaluated at complex argument U and
        %   parameter M=exp(-2*pi*L), 0 < L < Inf.  Recall that M = k^2,
        %   where k is the elliptic modulus.
        %
        %   U may be a matrix; L must be a scalar.  The entries of U are
        %   expected to lie within the rectangle |Re U| < K, 0 < Im U <
        %   Kp, where [K,Kp] = ELLIPK(L).
        %
        %   Copyright (c) 1999 by Toby Driscoll.
        %   $Id: ellipjc.m 298 2009-09-15 14:36:37Z driscoll $
        
        %   The built-in ELLIPJ can't handle compelx arguments, and
        %   standard transformations to handle this would require ELLIPJ
        %   called with parameter 1-M. When M < eps (or is even close),
        %   this can't be done accurately.
        %
        %   The algorithm is the descending Landen transformation,
        %   described in L. Howell's PhD thesis from MIT. Additional
        %   formulas from Gradshteyn & Ryzhik, 5th ed., and Abramowitz
        %   & Stegun.
        
        if nargin < 3
            % Absence of flag parameter indicates we must check for and transform u in
            % the upper half of the rectangle.
            [K,Kp] = rectmap.ellipkkp(L);
            high = imag(u) > Kp/2;
            u(high) = i*Kp - u(high);
            m = exp(-2*pi*L);
        else
            % Recursive call--L is actually m.
            high = zeros(size(u));
            m = L;
        end
        
        if m < 4*eps
            sinu = sin(u);
            cosu = cos(u);
            sn = sinu + m/4*(sinu.*cosu-u).*cosu;
            cn = cosu + m/4*(-sinu.*cosu+u).*sinu;
            dn = 1 + m/4*(cosu.^2-sinu.^2-1);
        else
            if m > 1e-3
                kappa = (1-sqrt(1-m))/(1+sqrt(1-m));
            else
                kappa = polyval([132,42,14,5,2,1,0],m/4);
            end
            mu = kappa^2;
            v = u/(1+kappa);
            [sn1,cn1,dn1] = rectmap.ellipjc(v,mu,1);
            denom = (1+kappa*sn1.^2);
            sn = (1+kappa)*sn1 ./ denom;
            cn = cn1.*dn1 ./ denom;
            dn = (1-kappa*sn1.^2) ./ denom;
        end
        
        if any(high(:))
            snh = sn(high);
            cnh = cn(high);
            dnh = dn(high);
            sn(high) = -1./(sqrt(m)*snh);
            cn(high) = i*dnh./(sqrt(m)*snh);
            dn(high) = i*cnh./snh;
        end
    end
    
    function [K,Kp] = ellipkkp(L)
        %ELLIPKKP Complete elliptic integral of the first kind, with complement.
        %   K = ELLIPKKP(L) returns the value of the complete elliptic
        %   integral of the first kind, evaluated at M=exp(-2*pi*L), 0 < L
        %   < Inf.
        %
        %   [K,KP] = ELLIPKKP(L) also returns the result for complementary
        %   parameter 1-M, which is useful when M < EPS.  Even when M <
        %   1e-6, the built-in ELLIPKE can lose digits of accuracy for KP.
        %
        %   Recall that the elliptic modulus k is related to the parameter
        %   M by M = k^2.
        %
        %   Copyright (c)1999 by Toby Driscoll.
        %   $Id: ellipkkp.m 298 2009-09-15 14:36:37Z driscoll $
        
        %   ELLIPKKP uses the method of the arithmetic-geometric mean described
        %   in 17.6 of M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
        %   Functions," Dover, 1965.  Same method as in ELLIPKE, only
        %   interchanging 1 and 1-m to find KP.
        
        % When m=exp(-2*pi*L) is extremely small, use O(m) approximations.
        if L > 10
            K = pi/2;
            Kp = pi*L + log(4);
            return
        end
        
        m = exp(-2*pi*L);
        a0 = 1;
        b0 = sqrt(1-m);
        s0 = m;
        i1 = 0; mm = 1;
        while mm > eps
            a1 = (a0+b0)/2;
            b1 = sqrt(a0.*b0);
            c1 = (a0-b0)/2;
            i1 = i1 + 1;
            w1 = 2^i1*c1.^2;
            mm = max(max(w1));
            s0 = s0 + w1;
            a0 = a1;
            b0 = b1;
        end
        K = pi./(2*a1);
        
        im = find(m==1);
        if ~isempty(im)
            K(im) = K(im)*inf;
        end
        
        if nargout > 1
            a0 = 1;
            b0 = sqrt(m);
            s0 = 1-m;
            i1 = 0; mm = 1;
            while mm > eps
                a1 = (a0+b0)/2;
                b1 = sqrt(a0.*b0);
                c1 = (a0-b0)/2;
                i1 = i1 + 1;
                w1 = 2^i1*c1.^2;
                mm = max(max(w1));
                s0 = s0 + w1;
                a0 = a1;
                b0 = b1;
            end
            Kp = pi./(2*a1);
            im = find(m==0);
            if ~isempty(im)
                Kp(im) = Kp(im)*inf;
            end
        end
    end
    
    function [yp,yprime] = r2strip(zp,z,L)
        %R2STRIP Map from rectangle to strip.
        %   R2STRIP(ZP,Z,L) maps from a rectangle to the strip 0 <= Im z <= 1,
        %   with the function log(sn(z|m))/pi, where sn is a Jacobi elliptic
        %   function and m = exp(-2*pi*L).  The prevertices of the map (in the
        %   rectangle domain) are given by Z; only the corners of the rectangle
        %   defined by Z are used.
        %
        %   The derivative of the map is returned as a second argument.
        %
        %   NOTE: The functionality is NOT parallel to HP2DISK and DISK2HP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: r2strip.m 298 2009-09-15 14:36:37Z driscoll $
        
        % The built-in ellipj accepts only real arguments. The standard identity is
        % not helpful when m is near zero, because (1-m) loses digits of accuracy.
        
        K = max(real(z));
        Kp = max(imag(z));
        yp = zp;
        yprime = zp;
        
        [sn,cn,dn] = diskmap.ellipjc(zp,L);
        % Make sure everything is in the upper half-plane (fix roundoff)
        sn = real(sn) + i*max(imag(sn),0);
        
        yp = log(sn)/pi;
        yprime = cn.*dn./sn/pi;
        
        % Make sure everything is in the strip (roundoff could put it outside)
        yp = real(yp) + i*max(0,imag(yp));
        yp = real(yp) + i*min(1,imag(yp));
    end
    
    function [w,beta,z,corners,renum] = rcorners(w,beta,z)
        %RCORNERS (not intended for calling directly by the user)
        %   Find corners of rectangle whose map is represented by prevertices z
        %   on the strip, then renumber w, beta, and z (and the corners) so that
        %   corners(1)=1.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: rcorners.m 298 2009-09-15 14:36:37Z driscoll $
        
        n = length(w);
        
        % Deduce corner locations
        left = abs(real(z)-min(real(z))) < eps;
        right = abs(real(z)-max(real(z))) < eps;
        top = abs(imag(z)-max(imag(z))) < eps;
        bot = abs(imag(z)-min(imag(z))) < eps;
        corners = find(left+right+top+bot - 1);
        c1 = find(abs(z-max(real(z))) < eps);
        offset = find(corners==c1);
        corners = corners([offset:4,1:offset-1]);
        
        % Renumber vertices so that corners(1)=1
        renum = [corners(1):n,1:corners(1)-1];
        w = w(renum);
        beta = beta(renum);
        z = z(renum);
        corners = rem(corners-corners(1)+1+n-1,n)+1;
    end
    
    function fprime = rderiv(zp,z,beta,c,L,zs)
        %RDERIV Derivative of the rectangle map.
        %   RDERIV(ZP,Z,BETA,C,L) returns the derivative at the points of ZP of
        %   the Schwarz-Christoffel rectangle map defined by Z, BETA, C, and L.
        %
        %   If a sixth argument is supplied, it is assumed to be the image of Z
        %   on the intermediate strip; see R2STRIP.
        %
        %   See also RPARAM, RMAP, R2STRIP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: rderiv.m 298 2009-09-15 14:36:37Z driscoll $
        
        n = length(z);
        
        if nargin < 6
            % Find prevertices on the strip
            zs = rectmap.r2strip(z,z,L);
            zs = real(zs) + i*round(imag(zs)); 	% put them *exactly* on edges
        end
        
        % First compute map and derivative from rectangle to strip
        [F,dF] = rectmap.r2strip(zp,z,L);
        
        % Now compute derivative of map from strip to polygon
        % Add in ends of strip
        ends = find(diff(imag(z([1:n 1]))));
        zs = [zs(1:ends(1));Inf;zs(ends(1)+1:ends(2));-Inf;zs(ends(2)+1:n)];
        bs = [beta(1:ends(1));0;beta(ends(1)+1:ends(2));0;beta(ends(2)+1:n)];
        dG = stripmap.deriv(F,zs,bs);
        
        % Put it together
        fprime = c*dF.*dG;
    end
    
    function rdisp(w,beta,z,c,L)
        %RDISP Display results of Schwarz-Christoffel rectangle parameter problem.
        %   RDISP(W,BETA,RECT,Z,C) displays the results of RPARAM in a pleasant
        %   way.
        %
        %   See also RPARAM, RPLOT.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: rdisp.m 298 2009-09-15 14:36:37Z driscoll $
        
        n = length(w);
        % Deduce corner locations
        left = abs(real(z)-min(real(z))) < eps;
        right = abs(real(z)-max(real(z))) < eps;
        top = abs(imag(z)-max(imag(z))) < eps;
        bot = abs(imag(z)-min(imag(z))) < eps;
        corners = find(left+right+top+bot - 1);
        c1 = find(abs(z-max(real(z))) < eps);
        offset = find(corners==c1);
        corners = corners([offset:4,1:offset-1]);
        rect = z(corners);
        
        disp(' ')
        disp(' cnr      vertex [w]          beta                prevertex [z]   ')
        disp(' ------------------------------------------------------------------------')
        u = real(w);
        v = imag(w);
        for j = 1:length(w)
            if v(j) < 0
                s = '-';
            else
                s = '+';
            end
            cnr = find(j==corners);
            if isempty(cnr)
                cstr = '    ';
            else
                cstr = sprintf('  %i ',cnr);
            end
            if ~imag(z(j))
                disp(sprintf('%s %8.5f %c %7.5fi    %8.5f   %16.8e',...
                    cstr,u(j),s,abs(v(j)),beta(j),z(j)));
            else
                disp(sprintf('%s %8.5f %c %7.5fi    %8.5f   %16.8e + %14.8ei',...
                    cstr,u(j),s,abs(v(j)),beta(j),real(z(j)),imag(z(j))));
            end
        end
        disp(' ')
        if imag(c) < 0
            s = '-';
        else
            s = '+';
        end
        disp(sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c))))
        disp(sprintf('\n  Conformal modulus = %.8g\n',imag(rect(2))/rect(1)/2));
    end
    
    function zp = rinvmap(wp,w,beta,z,c,L,qdat,z0,options)
        %RINVMAP Schwarz-Christoffel rectangle inverse map.
        %   RINVMAP(WP,W,BETA,CORNERS,Z,C,L,TOL) computes the inverse of the
        %   Schwarz-Christoffel rectangle map (i.e., from the polygon to the
        %   rectangle) at the points given in vector WP. The other arguments are
        %   as in RPARAM. TOL is a scalar tolerance, or a quadrature-data matrix
        %   QDAT as returned by SCQDATA, or may be omitted.
        %
        %   The default algorithm is to solve an ODE in order to obtain a fair
        %   approximation for ZP, and then improve ZP with Newton iterations.
        %   The ODE solution at WP requires a vector Z0 whose forward image W0
        %   is such that for each j, the line segment connecting WP(j) and W0(j)
        %   lies inside the polygon. By default Z0 is chosen by a fairly robust
        %   automatic process. Using a parameter (see below), you can choose to
        %   use either an ODE solution or Newton iterations exclusively.
        %
        %   RINVMAP(WP,...,TOL,Z0) has two interpretations.  If the ODE solution
        %   is being used, Z0 overrides the automatic selection of initial
        %   points. (This can be handy in convex polygons, where the choice of
        %   Z0 is trivial.) Otherwise, Z0 is taken as an initial guess to ZP. In
        %   either case, if length(Z0)==1, the value Z0 is used for all elements
        %   of WP; otherwise, length(Z0) should equal length(WP).
        %
        %   RINVMAP(WP,...,TOL,Z0,OPTIONS) uses a vector of parameters that
        %   control the algorithm. See SCINVOPT.
        %
        %   See also SCINVOPT, RPARAM, RMAP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: rinvmap.m 298 2009-09-15 14:36:37Z driscoll $
        
        import sctool.*
        n = length(w);
        w = w(:);
        beta = beta(:);
        z = z(:);
        [w,beta,z,corners] = rectmap.rcorners(w,beta,z);
        rect = z(corners);
        rect = [min(real(rect)) max(real(rect)) min(imag(rect)) max(imag(rect))];
        K = max(real(z));
        Kp = max(imag(z));
        zs = rectmap.r2strip(z,z,L);
        zs = real(zs) + 1i*round(imag(zs));	% put them *exactly* on edges
        
        zp = zeros(size(wp));
        wp = wp(:);
        lenwp = length(wp);
        
        if nargin < 9
            options = [];
            if nargin < 8
                z0 = [];
                if nargin < 7
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
                map = @(zp) rectmap.rmap(zp,w,beta,z,c,L,qdat);
                [z0,w0] = findz0('r',wp(~done),map,w,beta,z,c,L,qdat);
            else
                w0 = rectmap.rmap(z0,w,beta,z,c,L,qdat);
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
            odefun = @(w,y) rectmap.rimapfun(w,y,scale,z,beta,c,zs,L);
            [t,y] = ode23(odefun,[0,0.5,1],z0,odeset('abstol',odetol));
            [m,leny] = size(y);
            zp(~done) = y(m,1:lenwp)+sqrt(-1)*y(m,lenwp+1:leny);
            zp(~done) = rectmap.rectproject(zp(~done),rect);
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
                F = wp(~done) - rectmap.rmap(zn(~done),w,beta,z,c,L,qdat);
                dF = rectmap.rderiv(zn(~done),z,beta,c,L,zs);
                zn(~done) = zn(~done) + F(:)./dF(:);
                zn(~done) = rectmap.rectproject(zn(~done),rect);
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

    function zz = rectproject(zp,rect)
        % Project points into the rectangle.
        zz = max(min(real(zp),rect(2)),rect(1)) ...
            + 1i*max(min(imag(zp),rect(4)),rect(3));
    end

    function zdot = rimapfun(wp,yp,scale,z,beta,c,zs,L)
        %   Used by RINVMAP for solution of an ODE.

        %   Copyright 1998 by Toby Driscoll.
        %   $Id: rimapfun.m 298 2009-09-15 14:36:37Z driscoll $

        lenyp = length(yp);
        lenzp = lenyp/2;
        zp = yp(1:lenzp) + i*yp(lenzp+1:lenyp);

        f = scale./rectmap.rderiv(zp,z,beta,c,L,zs);
        zdot = [real(f);imag(f)];
    end
    
% Doesn't seem to be in use by the original SCT code, but had its own
% m-file. -- EK
%     function zdot = rimapfun(wp,yp,flag,scale,z,beta,c,zs,L);
%         %   Used by RINVMAP for solution of an ODE.
%         
%         %   Copyright 1998 by Toby Driscoll.
%         %   $Id: rimapfun.m 298 2009-09-15 14:36:37Z driscoll $
%         
%         lenyp = length(yp);
%         lenzp = lenyp/2;
%         zp = yp(1:lenzp) + i*yp(lenzp+1:lenyp);
%         
%         f = scale./rectmap.rderiv(zp,z,beta,c,L,zs);
%         zdot = [real(f);imag(f)];
%     end

    function wp = rmap(zp,w,beta,z,c,L,qdat)
        %RMAP   Schwarz-Christoffel rectangle map.
        %   RMAP(ZP,W,BETA,Z,C,L,QDAT) computes the values of the
        %   Schwarz-Christoffel rectangle map at the points in vector ZP.  The
        %   remaining arguments are as in RPARAM.  RMAP returns a vector the
        %   same size as ZP.
        %
        %   RMAP(ZP,W,BETA,Z,C,L,TOL) uses quadrature data intended to give an
        %   answer accurate to within TOL.
        %
        %   RMAP(ZP,W,BETA,Z,C,L) uses a tolerance of 1e-8.
        %
        %   See also RPARAM, RPLOT, RINVMAP.

        %   Copyright 1998 by Toby Driscoll.
        %   $Id: rmap.m 298 2009-09-15 14:36:37Z driscoll $

        if isempty(zp)
            wp = [];
            return
        end

        import sctool.*
        n = length(w);
        wp = z;
        w = w(:);
        beta = beta(:);
        z = z(:);
        [w,beta,z,corners] = rectmap.rcorners(w,beta,z);

        if nargin < 7
            qdat = scqdata(beta,8);
        elseif length(qdat)==1
            qdat = scqdata(beta,max(ceil(-log10(qdat)),8));
        end

        % Map prevertices to strip
        K = max(real(z));
        Kp = max(imag(z));
        zs = rectmap.r2strip(z,z,L);
        zs = real(zs) + i*round(imag(zs));	% put them *exactly* on edges

        % Add in ends of strip
        ends = find(diff(imag(zs([1:n 1]))));
        zs = [zs(1:ends(1));Inf;zs(ends(1)+1:ends(2));-Inf;zs(ends(2)+1:n)];
        ws = [w(1:ends(1));NaN;w(ends(1)+1:ends(2));NaN;w(ends(2)+1:n)];
        bs = [beta(1:ends(1));0;beta(ends(1)+1:ends(2));0;beta(ends(2)+1:n)];
        % Extend qdat with useless columns at ends
        idx = [1:ends(1) n+1 ends(1)+1:ends(2) n+1 ends(2)+1:n n+1];
        qdat = qdat(:,[idx idx+n+1]);

        wp = zeros(size(zp));
        zp = zp(:);
        p = length(zp);

        % Trap points which map to +/-Inf on the strip.
        bad = abs(zp) < 2*eps;
        zp(bad) = zp(bad) + 100*eps;
        bad = abs(zp-i*Kp) < 2*eps;
        zp(bad) = zp(bad) - i*100*eps*Kp;

        % Map from rectangle to strip.
        yp = rectmap.r2strip(zp,z,L);

        % Now map from strip to polygon.
        wp(:) = stripmap.evaluate(yp,ws,bs,zs,c,qdat);
    end
    
    function [z,c,L,qdat] = rparam(w,beta,cnr,z0,options)
        %RPARAM Schwarz-Christoffel rectangle parameter problem.
        %   [Z,C,L,QDAT] = RPARAM(W,BETA,CORNERS) solves the Schwarz-Christoffel
        %   parameter problem with a rectangle as fundamental domain and
        %   interior of the specified polygon as the target. W must be a vector
        %   of the vertices of the polygon, specified in counterclockwise
        %   order. BETA is a vector of turning angles; see SCANGLES. CORNERS is
        %   a 4-component vector specifying the indices of the vertices which
        %   are the images of the corners of the rectangle. **Be sure** the
        %   first two entries describe the LONG sides of the rectangle, and go
        %   in counterclockwise order. If CORNERS is omitted, the user is
        %   requested to select these vertices using the mouse.
        %
        %   If successful, RPARAM will return Z, a vector of the prevertices; C,
        %   the multiplicative constant of the conformal map; L, a parameter
        %   related to aspect ratio; and QDAT, an optional matrix of quadrature
        %   data used by some of the other SC routines.
        %
        %   [Z,C,L,QDAT] = RPARAM(W,BETA,CORNERS,Z0) uses Z0 as an initial guess
        %   for Z.  In this case, Z0 represents the image of prevertices on the
        %   strip 0 <= Im z <= 1.  You can use R2STRIP to transform prevertices
        %   from the rectangle to the strip.
        %
        %   [Z,C,L,QDAT] = RPARAM(W,BETA,CORNERS,TOL) attempts to find an answer
        %   within tolerance TOL. (Also see SCPAROPT.)
        %
        %   [Z,C,L,QDAT] = RPARAM(W,BETA,CORNERS,Z0,OPTIONS) uses a vector of
        %   control parameters. See SCPAROPT.
        %
        %   See also SCPAROPT, DRAWPOLY, RDISP, RPLOT, RMAP, RINVMAP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: rparam.m 298 2009-09-15 14:36:37Z driscoll $
        
        n = length(w); 				% no. of vertices
        w = w(:);
        beta = beta(:);
        import sctool.*
        
        % Set up defaults for missing args
        if nargin < 5
            options = [];
            if nargin < 4
                z0 = [];
                if nargin < 3
                    cnr = [];
                end
            end
        end
        
        [trace,tol,method] = parseopt(options);
        
        if isempty(cnr)
            msg{1} = 'Select the images of the corners of the rectangle.';
            msg{2} = 'Go in counterclockwise order and select a long rectangle edge first.';
            cnr = scselect(w,beta,4,'Select corners',msg);
        end
        
        % Renumber vertices so that cnr(1)=1
        renum = [cnr(1):n,1:cnr(1)-1];
        w = w(renum);
        beta = beta(renum);
        cnr = rem(cnr-cnr(1)+1+n-1,n)+1;
        
        if length(z0)==1
            tol = z0;
            z0 = [];
        end
        nqpts = max(ceil(-log10(tol)),4);
        qdat = scqdata(beta,nqpts); 		% quadrature data
        
        % Check input data.
        err = sccheck('r',w,beta,cnr);
        if err==1
            fprintf('Use SCFIX to make polygon obey requirements\n')
            error(' ')
        end
        
        atinf = (beta <= -1);
        
        if isempty(z0)
            % Try to find a reasonable initial guess.
            dw = abs(diff(w([1:n,1])));		% side lengths
            dw(isinf(dw)) = mean(dw(~isinf(dw)))*ones(size(dw(isinf(dw))));
            % Estimate length and width (conformal modulus)
            len = mean([sum(dw(cnr(1):cnr(2)-1)), sum(dw(cnr(3):cnr(4)-1))]);
            wid = mean([sum(dw(cnr(2):cnr(3)-1)), sum(dw([cnr(4):n,1:cnr(1)-1]))]);
            modest = min(len/wid,100);
            % Evenly space prevertices to match this conformal modulus
            z0(cnr(1):cnr(2)) = linspace(0,modest,diff(cnr(1:2))+1);
            dx = z0(cnr(1)+1)-z0(cnr(1));
            z0(cnr(1)-1:-1:1) = z0(cnr(1))-dx*(1:cnr(1)-1);
            z0(cnr(2)+1:cnr(3)-1) = z0(cnr(2)) + dx*(1:diff(cnr(2:3))-1);
            z0(cnr(4):-1:cnr(3)) = 1i + linspace(0,modest,diff(cnr(3:4))+1);
            dx = z0(cnr(4)-1)-z0(cnr(4));
            z0(cnr(4)+1:n) = z0(cnr(4))-dx*(1:n-cnr(4));
            
        else
            if length(z0)~=n
                error('Initial guess has wrong number of prevertices')
            end
            z0 = z0(renum);
            z0 = real(z0) + 1i*round(imag(z0));
            if any(imag(z0(cnr(1):cnr(2)))) || any(imag(z0(cnr(3):cnr(4)))==0)
                error('Initial guess has prevertices on wrong side of strip')
            end
        end
        
        % Convert z0 to unconstrained vars
        y0 = zeros(n-3,1);
        dz = diff(z0);
        dz(cnr(3):n-1) = -dz(cnr(3):n-1);
        y0(1:cnr(2)-2) = log(dz(1:cnr(2)-2));
        y0(cnr(3)-1:cnr(4)-3) = log(dz(cnr(3)+1:cnr(4)-1));
        y0(cnr(2)-1) = mean(log(dz([cnr(2)-1,cnr(3)])));
        
        % Vertices on the "short" edges are transformed into the interval [-1,1],
        % and then the Trefethen-style transformation is used.
        L = z0(cnr(2)) - z0(cnr(1));
        x = real(exp(pi*(L-conj(z0(cnr(2)+1:cnr(3)-1)))));
        dx = -diff([1;x(:);-1]);
        y0(cnr(2):cnr(3)-2) = log(dx(1:end-1)./dx(2:end));
        x = real(exp(pi*z0(cnr(4)+1:n)));
        dx = diff([-1;x(:);1]);
        y0(cnr(4)-2:n-3) = log(dx(1:end-1)./dx(2:end));
        
        
        % Find prevertices (solve param problem)
        
        % Set up normalized lengths for nonlinear equations:
        % indices of left and right integration endpoints
        left = 1:n-2;
        right = 2:n-1;
        % delete indices corresponding to vertices at Inf
        left(atinf(left)) = [];
        right(atinf(right)) = [];
        if atinf(n-1)
            right = [right,n];
        end
        cmplx = ((right-left) == 2);
        % It's possible we replaced the last single condition by a complex one.
        if length(cmplx)+sum(cmplx) > n-2
            cmplx(end) = 0;
        end
        % normalize lengths by w(2)-w(1)
        nmlen = (w(right)-w(left))/(w(2)-w(1));
        % abs value for finite ones
        nmlen(~cmplx) = abs(nmlen(~cmplx));
        % first entry is useless (=1)
        nmlen(1) = [];
        
        % Solve nonlinear system of equations:
        % package data
        %%figure
        %%hp = line(NaN,NaN);
        hp = [];
        fdat = {n,beta,nmlen,left,right,logical(cmplx),qdat,cnr};
        % set options
        opt = zeros(16,1);
        opt(1) = trace;
        opt(2) = method;
        opt(6) = 100*(n-3);			% max # of iterations
        opt(8) = tol;
        opt(9) = min(eps^(2/3),tol/10);
        opt(12) = nqpts;
        try
            [y,termcode] = nesolve(@rectmap.rpfun,y0,opt,fdat);
        catch
            % Have to delete the "waitbar" figure if interrupted
            close(findobj(allchild(0),'flat','Tag','TMWWaitbar'));
            error(lasterr)
        end
        if termcode~=1
            warning('Nonlinear equations solver did not terminate normally.')
        end
        
        % Convert y values to z on strip
        z = rectmap.rptrnsfm(y,cnr);
        ends = find(diff(imag(z([1:n 1]))));
        zs = [z(1:ends(1));Inf;z(ends(1)+1:ends(2));-Inf;z(ends(2)+1:n)];
        bs = [beta(1:ends(1));0;beta(ends(1)+1:ends(2));0;beta(ends(2)+1:n)];
        idx = [1:ends(1) n+1 ends(1)+1:ends(2) n+1 ends(2)+1:n n+1];
        qs = qdat(:,[idx idx+n+1]);
        
        % Determine multiplicative constant
        mid = mean(zs(1:2));
        g = stripmap.quad(zs(2),mid,2,zs,bs,qs) - stripmap.quad(zs(1),mid,1,zs,bs,qs);
        c = -diff(w(1:2))/g;
        
        
        % Find prevertices on the rectangle
        
        % Find corners of rectangle
        zs = z;
        L = max( real(zs(cnr([2 3]))) ) - zs(cnr(1));
        [K,Kp] = rectmap.ellipkkp(L);
        rect = [K;K+1i*Kp;-K+1i*Kp;-K];
        
        l = false(n,1);		% on left side
        l(cnr(3):cnr(4)) = 1;
        r = false(n,1);		% on right side
        r(cnr(1):cnr(2)) = 1;
        tl = (real(zs) > L) & imag(zs)>0;	% top-left side
        tr = (real(zs) > L) & imag(zs)==0;	% top-right side
        bl = (real(zs) < 0) & imag(zs)>0;	% bottom-left side
        br = (real(zs) < 0) & imag(zs)==0;	% bottom-right side
        
        % Initial guesses
        z = zeros(size(zs));
        % Corners
        z(cnr) = rect;
        % Left and right sides are simple
        z(l) = linspace(rect(3),rect(4),diff(cnr(3:4))+1).';
        z(r) = linspace(rect(1),rect(2),diff(cnr(1:2))+1).';
        % Cluster the top and bottom guesses near 0
        h = K/20;
        z(tl) = 1i*Kp - h*(1:sum(tl));
        z(tr) = 1i*Kp + h*(sum(tr):-1:1);
        z(bl) = -h*(sum(bl):-1:1);
        z(br) = h*(1:sum(br));
        
        % Newton iteration on "r2strip() - zs = 0"
        zn = z(:);
        maxiter = 50;
        done = false(size(zn));
        done(cnr) = 1;
        k = 0;  F = 0;
        while ~all(done) && k < maxiter
            [F,dF] = rectmap.r2strip(zn(~done),z(cnr),L);
            F = zs(~done) - F;
            
            % Adjust Newton step to stay exactly on rectangle boundary
            step = F./dF;
            lr = r(~done) | l(~done);			% on top or bottom
            step(lr) = 1i*imag(step(lr));
            step(~lr) = real(step(~lr));
            
            % Newton step
            znew = zn(~done) + step;
            
            % Keep prevertices from moving too far (past boundaries)
            % Left/right sides capped in Im direction
            x = min( max( imag(znew(lr)), 0 ), Kp);
            znew(lr) = real(znew(lr)) + i*x;
            % Top/bottom-left sides capped in Re direction
            tbl = tl(~done) | bl(~done);
            x = min( max(real(znew(tbl)), -K), -eps );
            znew(tbl) = i*imag(znew(tbl)) + x;
            % Top/bottom-right sides capped in Re direction
            tbr = tr(~done) | br(~done);
            x = min( max(real(znew(tbr)), eps), K );
            znew(tbr) = i*imag(znew(tbr)) + x;
            
            % Update
            zn(~done) = znew;
            done(~done) =  (abs(F) < tol);
            k = k + 1;
        end
        if any(abs(F)> tol)
            warning('Could not converge to the rectangle prevertices.')
        end
        z(:) = zn;
        
        % Undo renumbering
        z(renum) = z;
    end

    function F = rpfun(y,fdat)
        %   Returns residual for solution of nonlinear equations.
        %   $Id: rpfun.m 251 2003-03-07 16:34:11Z driscoll $
        
        [n,beta,nmlen,left,right,cmplx,qdat,corners] = deal(fdat{:});
        
        % Transform y (unconstr. vars) to z (actual params)
        z = rectmap.rptrnsfm(y,corners);
        
        % Compute the integrals appearing in nonlinear eqns.
        zleft = z(left);
        zright = z(right);
        
        % To use stquad, must put strip ends into z, beta
        ends = find(diff(imag(z([1:n 1]))));
        z = [z(1:ends(1));Inf;z(ends(1)+1:ends(2));-Inf;z(ends(2)+1:n)];
        beta = [beta(1:ends(1));0;beta(ends(1)+1:ends(2));0;beta(ends(2)+1:n)];
        % Put dummy columns into qdat at ends
        idx = [1:ends(1) n+1 ends(1)+1:ends(2) n+1 ends(2)+1:n n+1];
        qdat = qdat(:,[idx idx+n+1]);
        % Change singularity indices to reflect ends
        left = left + (left > ends(1)) + (left > ends(2));
        right = right + (right > ends(1)) + (right > ends(2));
        
        ints = zeros(size(zleft));
        % Two-stage integrations
        s2 = (right(:) - left(:) == 1) & (imag(zleft) - imag(zright) == 0);
        mid = mean([zleft(s2) zright(s2)],2);
        ints(s2) = stripmap.quadh(zleft(s2),mid,left(s2),z,beta,qdat) ...
            - stripmap.quadh(zright(s2),mid,right(s2),z,beta,qdat);
        
        % Three-stage integrations
        mid1 = real(zleft(~s2)) + 1i/2;
        mid2 = real(zright(~s2)) + 1i/2;
        ints(~s2) = stripmap.quad(zleft(~s2),mid1,left(~s2),z,beta,qdat) ...
            + stripmap.quadh(mid1,mid2,zeros(size(mid1)),z,beta,qdat) ...
            - stripmap.quad(zright(~s2),mid2,right(~s2),z,beta,qdat);
        
        if any(ints==0)|any(isnan(ints))
            % Singularities were too crowded.
            warning('Severe crowding')
        end
        
        % Compute nonlinear equation residual values.
        cmplx2 = cmplx(2:length(cmplx));
        F = abs(ints(~cmplx)); 		% F(1) = abs(ints(1))
        F = log( (F(2:end)/F(1)) ./ nmlen(~cmplx2) );
        if any(cmplx)
            F2 = log( (ints(cmplx)/ints(1)) ./ nmlen(cmplx2) );
            F = [F;real(F2);imag(F2)];
        end
    end
    
% Doesn't seem to be in use by the original SCT code, but it had its own
% m-file. -- EK
%     function F = rpfun(y,fdat)
%         %   Returns residual for solution of nonlinear equations.
%         %   $Id: rpfun.m 251 2003-03-07 16:34:11Z driscoll $
%         
%         [n,beta,nmlen,left,right,cmplx,qdat,corners] = deal(fdat{:});
%         
%         % Transform y (unconstr. vars) to z (actual params)
%         z = rectmap.rptrnsfm(y,corners);
%         
%         % Compute the integrals appearing in nonlinear eqns.
%         zleft = z(left);
%         zright = z(right);
%         
%         % To use stquad, must put strip ends into z, beta
%         ends = find(diff(imag(z([1:n 1]))));
%         z = [z(1:ends(1));Inf;z(ends(1)+1:ends(2));-Inf;z(ends(2)+1:n)];
%         beta = [beta(1:ends(1));0;beta(ends(1)+1:ends(2));0;beta(ends(2)+1:n)];
%         % Put dummy columns into qdat at ends
%         idx = [1:ends(1) n+1 ends(1)+1:ends(2) n+1 ends(2)+1:n n+1];
%         qdat = qdat(:,[idx idx+n+1]);
%         % Change singularity indices to reflect ends
%         left = left + (left > ends(1)) + (left > ends(2));
%         right = right + (right > ends(1)) + (right > ends(2));
%         
%         ints = zeros(size(zleft));
%         % Two-stage integrations
%         s2 = (right(:) - left(:) == 1) & (imag(zleft) - imag(zright) == 0);
%         mid = mean([zleft(s2) zright(s2)],2);
%         ints(s2) = stquadh(zleft(s2),mid,left(s2),z,beta,qdat) ...
%             - stquadh(zright(s2),mid,right(s2),z,beta,qdat);
%         
%         % Three-stage integrations
%         mid1 = real(zleft(~s2)) + i/2;
%         mid2 = real(zright(~s2)) + i/2;
%         ints(~s2) = stquad(zleft(~s2),mid1,left(~s2),z,beta,qdat) ...
%             + stquadh(mid1,mid2,zeros(size(mid1)),z,beta,qdat) ...
%             - stquad(zright(~s2),mid2,right(~s2),z,beta,qdat);
%         
%         if any(ints==0)|any(isnan(ints))
%             % Singularities were too crowded.
%             warning('Severe crowding')
%         end
%         
%         % Compute nonlinear equation residual values.
%         cmplx2 = cmplx(2:length(cmplx));
%         F = abs(ints(~cmplx)); 		% F(1) = abs(ints(1))
%         F = log( (F(2:end)/F(1)) ./ nmlen(~cmplx2) );
%         if any(cmplx)
%             F2 = log( (ints(cmplx)/ints(1)) ./ nmlen(cmplx2) );
%             F = [F;real(F2);imag(F2)];
%         end
%     end
    
    function [H,RE,IM] = rplot(w,beta,z,c,L,re,im,options)
        %RPLOT  Image of cartesian grid under Schwarz-Christoffel rectangle map.
        %   RPLOT(W,BETA,Z,C,L) will adaptively plot the images under the
        %   Schwarz-Christoffel rectangle map of ten evenly spaced horizontal
        %   and vertical lines in the retangle RECT.  The arguments are as in
        %   RPARAM.
        %
        %   RPLOT(W,BETA,Z,C,L,M,N) will plot images of M evenly spaced vertical
        %   and N evenly spaced horizontal lines.
        %
        %   RPLOT(W,BETA,Z,C,L,RE,IM) will plot images of vertical lines whose
        %   real parts are given in RE and horizontal lines whose imaginary
        %   parts are given in IM.  Either argument may be empty.
        %
        %   RPLOT(W,BETA,Z,C,L,RE,IM,OPTIONS) allows customization of RPLOT's
        %   behavior.  See SCPLTOPT.
        %
        %   H = RPLOT(W,BETA,Z,C,L,...) returns a vector of handles to all the
        %   curves drawn in the interior of the polygon.  [H,RE,IM] =
        %   RPLOT(W,BETA,Z,C,L,...) also returns the abscissae and ordinates of
        %   the lines comprising the grid.
        %
        %   See also SCPLTOPT, RPARAM, RMAP, RDISP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: rplot.m 298 2009-09-15 14:36:37Z driscoll $
        
        w = w(:);
        beta = beta(:);
        z = z(:);
        [w,beta,z,corners] = rectmap.rcorners(w,beta,z);
        rect = z(corners);
        import sctool.*
        
        % Parse input
        if nargin < 8
            options = [];
            if nargin < 7
                im = [];
                if nargin < 7
                    re = [];
                end
            end
        end
        
        Kp = imag(rect(2));
        K = rect(1);
        
        % Empty arguments default to 10
        if isempty([re(:);im(:)])
            re = 10;
            im = 10;
        end
        
        % Integer arguments must be converted to specific values
        if (length(re)==1) && (re == round(re))
            if re < 1
                re = [];
            else
                m = re;
                re = linspace(-K,K,m+2);
                re([1,m+2]) = [];
            end
        end
        if (length(im)==1) && (im == round(im))
            if im < 1
                im = [];
            else
                m = im;
                im = linspace(0,Kp,m+2);
                im([1,m+2]) = [];
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
        qdat = scqdata(beta,nqpts);
        len = max(diff(get(ax(1),'xlim')),diff(get(ax(1),'ylim')));
        minlen = len*minlen;
        maxlen = len*maxlen;
        axlim = axis;
        
        color = 'k';
        
        % Plot vertical lines...
        linh = gobjects(length(re),2);
        for j = 1:length(re)
            % Start evenly spaced
            zp = re(j) + 1i*linspace(0,Kp,15).';
            new = true(size(zp));
            wp = NaN(length(zp),1);
            
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
                neww = rectmap.rmap(zp(new),w,beta,z,c,L,qdat);
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
                    % Replace the points with the endpoints
                    set(linh(j,2),'erasemode','back')
                    set(linh(j,2),'marker','none','linestyle','-',...
                        'xdata',re(j)*[1 1],'ydata',[0 Kp])
                end
            else
                clearpoints(linh(j,1))
                addpoints(linh(j,1),real(wp(~isnan(wp))),imag(wp(~isnan(wp))));
                set(linh(j,1),'marker','none','linestyle','-','user',zp)
                if draw2
                    % Replace the points with smooth curve
                    clearpoints(linh(j,2))
                    addpoints(linh(j,2),re(j)*[1 1],[0 Kp])
                    set(linh(j,2),'marker','none','linestyle','-')
                end
            end
            drawnow
        end
        
        % Plot horizontal lines...
        linh1 = linh;
        linh = gobjects(length(im),2);
        for j = 1:length(im)
            % Start evenly spaced
            zp = linspace(-K,K,15).' + 1i*im(j);
            new = true(size(zp));
            wp = NaN(length(zp),1);
            
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
                neww = rectmap.rmap(zp(new),w,beta,z,c,L,qdat);
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
                    % Replace the points with the endpoints
                    set(linh(j,2),'erasemode','back')
                    set(linh(j,2),'marker','none','linestyle','-',...
                        'xdata',[-K K],'ydata',im(j)*[1 1])
                end
            else
                clearpoints(linh(j,1))
                addpoints(linh(j,1),real(wp(~isnan(wp))),imag(wp(~isnan(wp))));
                set(linh(j,1),'marker','none','linestyle','-','user',zp)
                if draw2
                    % Replace the points with (hopefully) a smooth circle
                    clearpoints(linh(j,2))
                    addpoints(linh(j,2),[-K K],im(j)*[1 1])
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
    
    function z = rptrnsfm(y,cnr)
        %RPTRNSFM (not intended for calling directly by the user)
        %   Transform optimization vars to prevertices for rectangle parameter
        %   problem.
        
        %       Copyright 1997 by Toby Driscoll. Last updated 05/06/97.
        
        n = length(y)+3;
        z = zeros(n,1);
        
        % Fill interior of "long edges" first
        z(cnr(1)+1:cnr(2)-1) = cumsum(exp(y(cnr(1):cnr(2)-2)));
        z(cnr(4)-1:-1:cnr(3)+1) = i + cumsum(exp(y(cnr(4)-3:-1:cnr(3)-1)));
        
        % Find L
        xr = real( z([cnr(2)-1,cnr(3)+1]) );
        z(cnr(2)) = mean(xr)+sqrt(diff(xr/2)^2+exp(2*y(cnr(2)-1)));
        z(cnr(3)) = i + z(cnr(2));
        z(cnr(4)) = i;
        
        % Now, fill in "short edges"
        cp = cumprod([1;exp(-y(cnr(2):cnr(3)-2))]);
        x = [0;cumsum(cp)] - [flipud(cumsum(flipud(cp)));0];
        x = x(2:end-1)/x(end);
        mask = abs(x) < eps;
        u = x;
        u(~mask) = log( x(~mask) ) / pi;
        u(mask) = -z(cnr(2))/eps;
        z(cnr(2)+1:cnr(3)-1) = i*imag(u) + real(z(cnr(2))) - real(u);
        
        idx = [cnr(4)-2:n-3 1:cnr(1)-1];
        cp = cumprod([1;exp(-y(idx))]);
        x = [0;cumsum(cp)] - [flipud(cumsum(flipud(cp)));0];
        x = x(2:end-1)/x(end);
        mask = abs(x) < eps;
        u = x;
        u(~mask) = log( x(~mask) ) / pi;
        u(mask) = -z(cnr(2))/eps;
        z([cnr(4)+1:n 1:cnr(1)-1]) = u;
    end
end

end
