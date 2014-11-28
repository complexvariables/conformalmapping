classdef stripmap < scmap
%STRIPMAP Schwarz-Christoffel strip map object.
%   STRIPMAP(P,ENDIDX) constructs a Schwarz-Christoffel strip map object
%   for the polygon P. ENDIDX is a two-vector containing the indices of
%   the vertices that are the images of the left and right ends of the
%   strip. STRIPMAP(P) requires you to choose these vertices
%   graphically. The parameter problem is solved using default options
%   for the prevertices and the multiplicative constant.
%
%   STRIPMAP(P,ENDIDX,OPTIONS) or STRIPMAP(P,OPTIONS) uses an options
%   structure of the type created by SCMAPOPT in solving the parameter
%   problem.
%
%   STRIPMAP(P,Z) creates a stripmap object having the given prevertices
%   Z (the mulitiplicative constant is found automatically). Z must
%   include -Inf and Inf, the ends of the strip.  STRIPMAP(P,Z,C) also
%   uses the given constant. An OPTIONS argument can be added, although
%   only the error tolerance will be used.
%
%   STRIPMAP(M), where M is a stripmap object, just returns M.
%
%   STRIPMAP(M,P) returns a new stripmap object for the polygon P using
%   the options in stripmap M. The prevertices of M will be used as the
%   starting guess for the parameter problem of the new map. Thus P
%   should properly be a perturbation (continuation) of the polygon for
%   M. An OPTIONS structure may be added to override options in M. There
%   is no opportunity to change the end indices.
%
%   See also SCMAPOPT, and classes POLYGON, SCMAP.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

%   Copyright 1998 by Toby Driscoll.

properties
    prevertex
    constant
    qdata
    accuracyVal
end

methods
    function map = stripmap(poly, varargin)
        % Assign empties to optional args
        endidx = [];
        z = [];
        c = [];
        opt = [];
        import sctool.*
        
        % Branch based on class of first argument
        switch class(poly)
            
            case 'stripmap'
                oldmap = poly;
                % Continuation of given map to given polygon
                poly = varargin{1};
                z0 = oldmap.prevertex;
                if length(z0) ~= length(poly)
                    msg = 'Polygon %s must have the same length as that in %s.';
                    error(msg,inputname(2),inputname(1))
                end
                opt = scmapopt(oldmap.scmap);
                if nargin > 2
                    opt = scmapopt(opt,varargin{2});
                end
                opt = scmapopt(opt,'initial',z0);
                endidx(1) = find( isinf(z0) & (z0 < 0) );
                endidx(2) = find( isinf(z0) & (z0 > 0) );
                
            case 'polygon'
                % Parse optional arguments
                for j = 1:length(varargin)
                    arg = varargin{j};
                    % Each arg is the end specifier, an options struct, z, or c
                    if isa(arg,'struct')
                        opt = arg;
                    elseif (length(arg) == 2) & all(round(arg) == arg)
                        endidx = arg;
                    elseif length(arg) == length(poly)
                        z = arg;
                        z = z(:);
                    elseif length(arg) == 1
                        c = arg;
                    else
                        msg = 'Unable to parse argument ''%s''.';
                        error(msg,inputname(j+1))
                    end
                end
                
            otherwise
                msg = 'Expected ''%s'' to be of class polygon or stripmap.';
                error(msg,inputname(1))
                
        end % switch
        
        % Retrieve options
        opt = scmapopt(opt);
        
        
        % Get data for the low-level functions
        w = vertex(poly);
        n = length(w);
        beta = angle(poly) - 1;
        
        % Request endidx
        if isempty(endidx)
            msg = 'Select the vertices that map to the ends of the strip.';
            endidx = scselect(w,beta,2,'Select ends',msg);
        end
        
        % Find prevertices if necessary
        if isempty(z)
            
            % Apply SCFIX to enforce solver rules
            [w,beta,endidx] = scfix('st',w,beta,endidx);
            poly = polygon(w,beta+1);
            
            [z,c,qdata] = stripmap.stparam(w,beta,endidx,opt.InitialGuess,opt);
            
        else
            
            atinf = isinf(z);
            if (sum(atinf & (z < 0)) ~= 1) | (sum(atinf & (z > 0)) ~= 1)
                error('Supplied prevertices must include -Inf and Inf')
            end
            % Base quadrature accuracy on given options
            nqpts = ceil(-log10(opt.Tolerance));
            qdata = scqdata(beta,nqpts);
            
        end
        
        % Find constant if necessary
        if isempty(c)
            k = find(isinf(z) & (z < 0));
            mid = mean(z(k+1:k+2));
            I = stripmap.stquad(z(k+1),mid,k+1,z,beta,qdata) - ...
                stripmap.stquad(z(k+2),mid,k+2,z,beta,qdata);
            c = diff(w(k+1:k+2))/I;
        end
        
        map = map@scmap(poly,opt);
        map.prevertex = z;
        map.constant = c;
        map.qdata = qdata;
        
        % Now fill in apparent accuracy
        map.accuracyVal = accuracy(map);
    end
    
    function acc = accuracy(M)
        %ACCURACY Apparent accuracy of Schwarz-Christoffel strip map.
        %   ACCURACY(M) estimates the accuracy of the Schwarz-Christoffel strip
        %   map M. The technique used is to compare the differences between
        %   successive finite vertices to the integral between the corresponding
        %   prevertices, and return the maximum.
        %
        %   See also STRIPMAP.
        
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
        qdata = M.qdata;
        
        % Find ends of strip and renumber to put left end first
        endidx(1) = find( isinf(z) & (z < 0) );
        endidx(2) = find( isinf(z) & (z > 0) );
        
        % Integrate between consecutive finite pairs along bottom and top
        bot = ~imag(z);
        bot(endidx) = 0;
        top = logical(imag(z));
        top(endidx) = 0;
        idxbot = find( bot & ~isinf(w) );
        idxtop = find( top & ~isinf(w) );
        
        % Two columns hold endpoint indices for integrations
        idx = [idxbot(1:end-1) idxbot(2:end)];
        idx = [idx ; [idxtop(1:end-1) idxtop(2:end)] ];
        
        % As a final check, integrate once across the strip
        [tmp,k] = min(abs( real(z(idxtop)-z(rem(endidx(1),n)+1)) ));
        idx = [idx; [rem(endidx(1),n)+1 idxtop(k)]];
        
        I = zeros(size(idx,1),1);
        zl = z(idx(:,1));
        zr = z(idx(:,2));
        
        % Two-stage integrations (neighboring prevertices)
        s2 = (diff(idx,1,2) == 1);
        mid = mean([zl(s2) zr(s2)],2);
        I(s2) = M.stquadh(zl(s2),mid,idx(s2,1),z,beta,qdata) ...
            - M.stquadh(zr(s2),mid,idx(s2,2),z,beta,qdata);
        
        % Three-stage integrations
        mid1 = real(zl(~s2)) + i/2;
        mid2 = real(zr(~s2)) + i/2;
        I(~s2) = M.stquad(zl(~s2),mid1,idx(~s2,1),z,beta,qdata) ...
            + M.stquadh(mid1,mid2,zeros(size(mid1)),z,beta,qdata) ...
            - M.stquad(zr(~s2),mid2,idx(~s2,2),z,beta,qdata);
        
        acc = max(abs( c*I - diff(w(idx).').' ));
    end
    
    function out = char(f)
        %CHAR   Pretty-print a Schwarz-Christoffel strip map.
        
        %   Copyright 1998-2001 by Toby Driscoll.
        %   $Id: char.m 162 2001-07-20 14:33:00Z driscoll $
        
        p = f.polygon;
        w = vertex(p);
        alpha = angle(p);
        z = f.prevertex;
        c = f.constant;
        
        L = cell(2,1);
        L{1} = '      vertex               alpha           prevertex         ';
        L{2}=' ------------------------------------------------------------';
        u = real(w);
        v = imag(w);
        for j = 1:length(w)
            if v(j) < 0
                s = '-';
            else
                s = '+';
            end
            if ~imag(z(j))
                L{end+1}=sprintf(' %8.5f %c %7.5fi     %8.5f    %20.12e',...
                    u(j),s,abs(v(j)),alpha(j),z(j));
            else
                L{end+1}=sprintf(' %8.5f %c %7.5fi     %8.5f    %20.12e + i',...
                    u(j),s,abs(v(j)),alpha(j),z(j));
            end
        end
        
        L{end+1} = ' ';
        if imag(c) < 0
            s = '-';
        else
            s = '+';
        end
        L{end+1} = sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c)));
        L{end+1} = sprintf('  Apparent accuracy is %.2e',f.accuracyVal);
        L{end+1} = ' ';
        
        out = L;
    end
    
    function Md = diskmap(Ms)
        %   Convert strip map to disk map.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: diskmap.m 43 1998-07-01 19:12:31Z tad $
        
        p = Ms.polygon;
        w = vertex(p);
        beta = angle(p) - 1;
        z = Ms.prevertex;
        n = length(z);
        
        % Find index of vertex at Inf
        idx = find(isinf(z) & (z > 0));
        
        % Put that vertex last
        renum = [idx+1:n 1:idx];
        w = w(renum);
        beta = beta(renum);
        z = z(renum);
        
        % Map prevertices to real axis
        zh = exp(pi*z);
        
        % Map -Inf correctly
        idx = find(isinf(z) & (z < 0));
        zh(idx) = 0;
        
        % Map Inf correctly
        zh(n) = Inf;
        
        % Transform prevertices to disk
        A = moebius(zh(n-2:n),[-1 -i 1]);
        zd = sign(A(zh));
        zd(n) = 1;
        
        % Create new map
        Md = diskmap(polygon(w,beta+1),zd);
    end
    
    function out = display(M)
        %DISPLAY Display parameters of a Schwarz-Christoffel strip map.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: display.m 118 2001-05-03 21:18:27Z driscoll $
        
        p = M.polygon;
        w = vertex(p);
        alpha = angle(p);
        z = M.prevertex;
        c = M.constant;
        
        L = {' '; '  stripmap object:'; ' '};
        L{4} = '      vertex               alpha           prevertex         ';
        L{5}=' ------------------------------------------------------------';
        u = real(w);
        v = imag(w);
        for j = 1:length(w)
            if v(j) < 0
                s = '-';
            else
                s = '+';
            end
            if ~imag(z(j))
                L{end+1}=sprintf(' %8.5f %c %7.5fi     %8.5f    %20.12e',...
                    u(j),s,abs(v(j)),alpha(j),z(j));
            else
                L{end+1}=sprintf(' %8.5f %c %7.5fi     %8.5f    %20.12e + i',...
                    u(j),s,abs(v(j)),alpha(j),z(j));
            end
        end
        
        L{end+1} = ' ';
        if imag(c) < 0
            s = '-';
        else
            s = '+';
        end
        L{end+1} = sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c)));
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
        %EVAL Evaluate Schwarz-Christoffel strip map at points.
        %   EVAL(M,ZP) evaluates the Schwarz-Christoffel map M at the points ZP
        %   in the strip 0<=Im z<=1. The default tolerance of M is used.
        %
        %   EVAL(M,ZP,TOL) attempts to give an answer accurate to TOL. If TOL is
        %   less than the accuracy of M, this is unlikely to be met.
        %
        %   See also STRIPMAP, EVALINV.
        
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
        idx = (imag(zp) > -eps) & (imag(zp) < 1+eps);
        wp(idx) = ...
            M.stmap(zp(idx),vertex(p),angle(p)-1,M.prevertex,M.constant,qdata);
    end
    
    function fp = evaldiff(M,zp)
        %EVALDIFF Derivative of Schwarz-Christoffel strip map at points.
        %   EVALDIFF(M,ZP) computes the derivative of the Schwarz-Christoffel
        %   strip map M at the points ZP.
        %
        %   See also STRIPMAP, EVAL.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: evaldiff.m 7 1998-05-10 04:37:19Z tad $
        
        z = M.prevertex;
        c = M.constant;
        beta = angle(M.polygon) - 1;
        
        fp = M.stderiv(zp,z,beta,c);
    end
    
    function zp = evalinv(M,wp,tol,z0)
        %EVALINV Invert Schwarz-Christoffel strip map at points.
        %   EVALINV(M,WP) evaluates the inverse of the Schwarz-Christoffel map M
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
        %   See also STRIPMAP, EVAL.
        
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
        
        zp = NaN*wp;
        %idx = logical(isinpoly(wp,p));
        idx = logical(ones(size(wp)));
        zp(idx) = M.stinvmap(wp(idx),w,beta,z,c,qdata,z0,[0 tol]);
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
    
    function Mh = hplmap(Ms)
        %   Convert strip map to half-plane map.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: hplmap.m 7 1998-05-10 04:37:19Z tad $
        
        p = Ms.polygon;
        w = vertex(p);
        beta = angle(p) - 1;
        z = Ms.prevertex;
        n = length(z);
        
        % Find index of vertex at Inf
        idx = find(isinf(z) & (z > 0));
        
        % Put that vertex last
        renum = [idx+1:n 1:idx];
        w = w(renum);
        beta = beta(renum);
        z = z(renum);
        
        % Map prevertices to real axis
        zh = real(exp(pi*z));
        
        % Map -Inf correctly
        idx = find(isinf(z) & (z < 0));
        zh(idx) = 0;
        
        % Put finite ones inside [-1,1]
        zh = (zh-zh(1)) * 2/(zh(n-1)-zh(1)) - 1;
        
        % Map Inf correctly
        zh(n) = Inf;
        
        % Create new map
        Mh = hplmap(polygon(w,beta+1),zh);
    end
    
    function M = mtimes(M,c)
        %   Scale the image of a map by a complex constant.
        
        %   Copyright (c) 1998 by Toby Driscoll.
        %   $Id: mtimes.m 33 1998-06-29 22:35:40Z tad $
        
        % May need to swap arguments
        if isa(M,'double') & isa(c,'stripmap')
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
        %PLOT Visualize a Schwarz-Christoffel strip map.
        %   PLOT(M) plots the polygon associated with the Schwarz-Christoffel
        %   strip map M and the images of ten evenly spaced vertical line
        %   segments and horizontal lines under the S-C transformation.
        %
        %   PLOT(M,NRE,NIM) plots the images of NRE vertical line segments and
        %   NIM horizontal lines.
        %
        %   PLOT(M,RE,IM) plots the vertical line segments at abscissae given by
        %   the entries of RE and horizontal lines at the ordinates specified in
        %   IM.
        %
        %   PLOT(M,TOL) or PLOT(M,NRE,NIM,TOL) or PLOT(M,RE,IM,TOL) computes the
        %   map with accuracy roughly TOL. Normally TOL defaults to 1e-4 or the
        %   accuracy of M, whichever is greater.
        %
        %   See also STRIPMAP, EVAL.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: plot.m 7 1998-05-10 04:37:19Z tad $
        
        p = M.polygon;
        w = vertex(p);
        beta = angle(p) - 1;
        z = M.prevertex;
        c = M.constant;
        
        if nargin == 1
            [a1,a2,a3] = M.stplot(w,beta,z,c);
        elseif length(varargin) == 1
            % Tolerance given only
            [a1,a2,a3] = M.stplot(w,beta,z,c,10,10,ceil(-log10(varargin{1})));
        elseif length(varargin) == 2
            % RE,IM given only
            [a1,a2,a3] = M.stplot(w,beta,z,c,varargin{1},varargin{2});
        else
            % All given
            nqpts = ceil(-log10(varargin{3}));
            [a1,a2,a3] = M.stplot(w,beta,z,c,varargin{1},varargin{2},nqpts);
        end
        
        if nargout > 0
            h = a1;
            re = a2;
            im = a3;
        end
    end
end

methods (Static)
    % Primitives for rectangle maps need these.
    
    function q = deriv(varargin)
        q = stripmap.stderiv(varargin{:});
    end
    
    function q = quad(varargin)
        q = stripmap.stquad(varargin{:});
    end
    
    function q = quadh(varargin)
        q = stripmap.stquadh(varargin{:});
    end
    
    function w = evaluate(varargin)
        w = stripmap.stmap(varargin{:});
    end
end

methods(Hidden,Static)
    function fprime = stderiv(zp,z,beta,c,j)
        %STDERIV Derivative of the strip map.
        %   STDERIV(ZP,Z,BETA,C) returns the derivative at the points of ZP of
        %   the Schwarz-Christoffel strip map defined by Z, BETA, and C.
        %
        %   See also STPARAM, STMAP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: stderiv.m 298 2009-09-15 14:36:37Z driscoll $
        
        %   If a fifth argument j is supplied, the terms corresponding to z(j)
        %   are normalized by abs(zp-z(j)).  This is for Gauss-Jacobi
        %   quadrature.
        
        % Support old syntax
        if nargin < 4
            c = 1;
        end
        
        log2 = 0.69314718055994531;
        fprime = zeros(size(zp));
        zprow = zp(:).';
        npts = length(zprow);
        
        % Strip out infinite prevertices
        if length(z)==length(beta)
            ends = find(isinf(z));
            theta = diff(beta(ends));
            if z(ends(1)) < 0
                theta = -theta;
            end
            z(ends) = [];
            beta(ends) = [];
            % Adjust singularity index if given
            if nargin > 4
                j = j - (j > ends(1)) - (j > ends(2));
            end
        else
            error('Vector of prevertices must include +/-Inf entries')
        end
        zcol = z(:);
        bcol = beta(:);
        n = length(z);
        
        terms = -pi/2*(zprow(ones(n,1),:) - zcol(:,ones(npts,1)));
        lower = (~imag(z));
        terms(lower,:) = -terms(lower,:);
        rt = real(terms);
        big = abs(rt) > 40;
        if any(any(~big))
            terms(~big) = log(-i*sinh(terms(~big)));
        end
        terms(big) = sign(rt(big)).*(terms(big)-i*pi/2) - log2;
        if nargin > 4
            if j > 0
                terms(j,:) = terms(j,:)-log(abs(zprow-z(j)));
            end
        end
        fprime(:) = c*exp(pi/2*theta*zprow + sum(terms.*bcol(:,ones(npts,1))));
    end
    
    function stdisp(w,beta,z,c)
        %STDISP Display results of Schwarz-Christoffel strip parameter problem.
        %   STDISP(W,BETA,Z,C) displays the results of STPARAM in a pleasant
        %   way.
        %
        %   See also STPARAM, STPLOT.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: stdisp.m 298 2009-09-15 14:36:37Z driscoll $
        
        disp(' ')
        disp('      vertex [w]           beta            prevertex [z]     ')
        disp(' ------------------------------------------------------------')
        u = real(w);
        v = imag(w);
        for j = 1:length(w)
            if v(j) < 0
                s = '-';
            else
                s = '+';
            end
            if ~imag(z(j))
                disp(sprintf(' %8.5f %c %7.5fi     %8.5f    %20.12e',...
                    u(j),s,abs(v(j)),beta(j),z(j)));
            else
                disp(sprintf(' %8.5f %c %7.5fi     %8.5f    %20.12e + i',...
                    u(j),s,abs(v(j)),beta(j),z(j)));
            end
        end
        disp(' ')
        if imag(c) < 0
            s = '-';
        else
            s = '+';
        end
        disp(sprintf('  c = %.8g %c %.8gi\n',real(c),s,abs(imag(c))))
    end
    
% Is this used? Had its own m-file. -- EK
%     function zdot = stimapfun(wp,yp,flag,scale,z,beta,c);
%         
%         %   Used by STINVMAP for solution of an ODE.
%         
%         %   Copyright 1998 by Toby Driscoll.
%         %   $Id: stimapfun.m 298 2009-09-15 14:36:37Z driscoll $
%         
%         lenyp = length(yp);
%         lenzp = lenyp/2;
%         zp = yp(1:lenzp) + i*yp(lenzp+1:lenyp);
%         
%         f = scale./stripmap.stderiv(zp,z,beta,c);
%         zdot = [real(f);imag(f)];
%     end
    
    function zp = stinvmap(wp,w,beta,z,c,qdat,z0,options)
        %STINVMAP Schwarz-Christoffel strip inverse map.
        %   STINVMAP(WP,W,BETA,Z,C,TOL) computes the inverse of the
        %   Schwarz-Christoffel strip map (i.e., from the polygon to the strip)
        %   at the points given in vector WP. The other arguments are as in
        %   STPARAM. TOL is a scalar tolerance, or a quadrature-data matrix QDAT
        %   as returned by SCQDATA, or may be omitted.
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
        %   STINVMAP(WP,W,BETA,Z,C,TOL,Z0) has two interpretations. If the ODE
        %   solution is being used, Z0 overrides the automatic selection of
        %   initial points. (This can be handy in convex polygons, where the
        %   choice of Z0 is trivial.) Otherwise, Z0 is taken as an initial guess
        %   to ZP. In either case, if length(Z0)==1, the value Z0 is used for
        %   all elements of WP; otherwise, length(Z0) should equal length(WP).
        %
        %   STINVMAP(WP,W,BETA,Z,C,TOL,Z0,OPTIONS) uses a vector of parameters
        %   that control the algorithm. See SCINVOPT.
        %
        %   See also SCINVOPT, STPARAM, STMAP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: stinvmap.m 298 2009-09-15 14:36:37Z driscoll $
        
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
                map = @(zp) stripmap.stmap(zp,w,beta,z,c,qdat);
                [z0,w0] = findz0('st',wp(~done),map,w,beta,z,c,qdat);
            else
                w0 = stripmap.stmap(z0,w,beta,z,c,qdat);
                if length(z0)==1 & lenwp > 1
                    z0 = z0(:,ones(lenwp,1)).';
                    w0 = w0(:,ones(lenwp,1)).';
                end
                w0 = w0(~done);
                z0 = z0(~done);
            end
            
            % Use relaxed ODE tol if improving with Newton.
            odetol = max(tol,1e-2*(newton));
            
            % Rescale dependent coordinate
            scale = (wp(~done) - w0(:));
            
            % Solve ODE
            z0 = [real(z0);imag(z0)];
            odefun = @(w,y) stripmap.stimapfun(w,y,scale,z,beta,c);
            [t,y] = ode23(odefun,[0,0.5,1],z0,odeset('abstol',odetol));
            [m,leny] = size(y);
            zp(~done) = y(m,1:lenwp)+sqrt(-1)*y(m,lenwp+1:leny);
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
                F = wp(~done) - stripmap.stmap(zn(~done),w,beta,z,c,qdat);
                dF = c*stripmap.stderiv(zn(~done),z,beta);
                zn(~done) = zn(~done) + F(:)./dF(:);
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

    function zdot = stimapfun(wp,yp,scale,z,beta,c)
        %   Used by STINVMAP for solution of an ODE.

        %   Copyright 1998 by Toby Driscoll.
        %   $Id: stimapfun.m 298 2009-09-15 14:36:37Z driscoll $

        lenyp = length(yp);
        lenzp = lenyp/2;
        zp = yp(1:lenzp) + i*yp(lenzp+1:lenyp);

        f = scale./stripmap.stderiv(zp,z,beta,c);
        zdot = [real(f);imag(f)];
    end
    
    function wp = stmap(zp,w,beta,z,c,qdat)
        %STMAP  Schwarz-Christoffel strip map.
        %   STMAP(ZP,W,BETA,Z,C,QDAT) computes the values of the Schwarz-
        %   Christoffel strip map at the points in vector ZP. The arguments W,
        %   BETA, Z, C, and QDAT are as in STPARAM. STMAP returns a vector the
        %   same size as ZP.
        %
        %   STMAP(ZP,W,BETA,Z,C,TOL) uses quadrature data intended to give an
        %   answer accurate to within TOL.
        %
        %   STMAP(ZP,W,BETA,Z,C) uses a tolerance of 1e-8.
        %
        %   See also STPARAM, STPLOT, STINVMAP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: stmap.m 298 2009-09-15 14:36:37Z driscoll $
        
        if isempty(zp)
            wp = [];
            return
        end
        
        import sctool.*
        n = length(w);
        w = w(:);
        beta = beta(:);
        z = z(:);
        
        % Quadrature data
        if nargin < 6
            qdat = scqdata(beta,8);
        elseif length(qdat)==1
            qdat = scqdata(beta,max(ceil(-log10(qdat)),8));
        end
        tol = 10^(-size(qdat,1));
        
        shape = size(zp);
        zp = zp(:);
        zprow = zp.';
        p = length(zp);
        wp = zeros(p,1);
        
        % For each point in zp, find the nearest prevertex.
        [dist,sing] = min( abs(zprow(ones(n,1),:) - z(:,ones(1,p))) );
        sing = sing(:);				% indices of prevertices
        
        % Screen out images of prevertices
        vertex = (dist(:) < tol);
        wp(vertex) = w(sing(vertex));
        leftend = (isinf(zp) & (real(zp) < 0));
        wp(leftend) = w(z == -Inf);
        rightend = (isinf(zp) & (real(zp) > 0));
        wp(rightend) = w(z == Inf);
        vertex = vertex | leftend | rightend;
        
        % "Bad" points are closest to a prevertex of infinity.
        atinf = find(isinf(w(:))); 		% infinite vertices
        bad = ismember(sing,atinf) & ~vertex;
        
        if any(bad)
            % We can't begin integrations at a preimage of infinity. We will pick the
            % next-closest qualified prevertex.
            zf = z;
            zf(isinf(w)) = Inf;
            [tmp,s] = min( abs(zprow(ones(n,1),bad) - zf(:,ones(sum(bad),1))) );
            sing(bad) = s;
            
            % Because we no longer integrate from the closest prevertex, we must go in
            % stages to maintain accuracy.
            mid1 = real(z(s)) + i/2;
            mid2 = real(zp(bad)) + i/2;
        else
            bad = zeros(p,1);		% all clear
        end
        
        % zs = the starting singularities
        zs = z(sing);
        % ws = map(zs)
        ws = w(sing);
        
        % Compute the map directly at "normal" points.
        normal = ~bad & ~vertex;
        if any(normal)
            I = stripmap.stquad(zs(normal),zp(normal),sing(normal),z,beta,qdat);
            wp(normal) = ws(normal) + c*I;
        end
        
        % Compute map at "bad" points in stages.
        if any(bad)
            I1 = stripmap.stquad(zs(bad),mid1,sing(bad),z,beta,qdat);
            I2 = stripmap.stquadh(mid1,mid2,zeros(sum(bad),1),z,beta,qdat);
            I3 = -stripmap.stquad(zp(bad),mid2,zeros(sum(bad),1),z,beta,qdat);
            wp(bad) = ws(bad) + c*(I1 + I2 + I3);
        end
        
        wp = reshape(wp,shape);
    end
    
    function [z,c,qdat] = stparam(w,beta,ends,z0,options)
        %STPARAM Schwarz-Christoffel strip parameter problem.
        %   [Z,C,QDAT] = STPARAM(W,BETA,ENDS) solves the Schwarz-Christoffel
        %   parameter problem with the infinite strip as fundamental domain and
        %   interior of the specified polygon as the target. W must be a vector
        %   of the vertices of the polygon, specified in counterclockwise
        %   order. BETA is a vector of turning angles; see SCANGLES. ENDS is a
        %   2-vector whose entries are the indices of the vertices which are the
        %   images of the left and right ends of the strip. If ENDS is omitted,
        %   the user is requested to select these vertices using the mouse.
        %
        %   If successful, STPARAM will return Z, a vector of the pre-images of
        %   W; C, the multiplicative constant of the conformal map; and QDAT, an
        %   optional matrix of quadrature data required by some of the other SC
        %   routines.
        %
        %   [Z,C,QDAT] = STPARAM(W,BETA,ENDS,Z0) uses Z0 as an initial guess for
        %   Z.
        %
        %   [Z,C,QDAT] = STPARAM(W,BETA,ENDS,TOL) attempts to find an answer
        %   within tolerance TOL. (Also see SCPAROPT.)
        %
        %   [Z,C,QDAT] = STPARAM(W,BETA,ENDS,Z0,OPTIONS) uses a vector of
        %   control parameters.  See SCPAROPT.
        %
        %   See also SCPAROPT, DRAWPOLY, STDISP, STPLOT, STMAP, STINVMAP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: stparam.m 298 2009-09-15 14:36:37Z driscoll $
        
        % Set up defaults for missing args
        if nargin < 5
            options = [];
            if nargin < 4
                z0 = [];
                if nargin < 3
                    ends = [];
                end
            end
        end
        import sctool.*
        [trace,tol,method] = parseopt(options);
        
        if isempty(ends)
            msg = 'Select the vertices that map to the ends of the strip.';
            ends = scselect(w,beta,2,'Select ends',msg);
        end
        
        N = length(w); 				% no. of vertices
        w = w(:);
        beta = beta(:);
        % Renumber vertices so that the ends of the strip map to w([1,k])
        renum = [ends(1):N,1:ends(1)-1];
        w = w(renum);
        beta = beta(renum);
        k = find(renum==ends(2));
        % n: number of finite prevertices
        n = N-2;
        % nb: number of prevertices on bottom edge of strip
        nb = k-2;
        
        % Check input data.
        err = sccheck('st',w,beta,ends);
        if err==1
            fprintf('Use SCFIX to make polygon obey requirements\n')
            error(' ')
        end
        
        if length(z0)==1
            tol = z0;
            z0 = [];
        end
        nqpts = max(ceil(-log10(tol)),4);
        qdat = scqdata(beta,nqpts); 		% quadrature data
        atinf = (beta <= -1);
        
        % Ignore images of ends of strip.
        w([1,k]) = [];
        atinf([1,k]) = [];
        
        if isempty(z0)
            % Make initial guess based on polygon.
            z0 = zeros(n,1);
            if any(atinf)
                % Can't base it on relative side lengths.
                scale = (abs(w(nb)-w(1))+abs(w(n)-w(nb+1)))/2;
                z0(1:nb) = linspace(0,scale,nb)';
                z0(nb+1:n) = i + flipud(linspace(0,scale,n-nb)');
            else
                % This is from Louis Howell's code.
                scale = (abs(w(n)-w(1))+abs(w(nb)-w(nb+1)))/2;
                z0(1:nb) = cumsum([0;abs(w(2:nb)-w(1:nb-1))]/scale);
                if nb+1==n
                    z0(n) = mean(z0([1,nb]));
                else
                    z0(n:-1:nb+1) = cumsum([0;abs(w(n:-1:nb+2)-w(n-1:-1:nb+1))]/scale);
                end
                scale = sqrt(z0(nb)/z0(nb+1));
                z0(1:nb) = z0(1:nb)/scale;
                z0(nb+1:n) = i + z0(nb+1:n)*scale;
            end
        else
            z0 = z0(renum);
            if length(z0)==N
                if ~all(isinf(z0([1,k])))
                    error('Starting guess does not match ends of strip')
                end
                z0([1,k]) = [];
            elseif length(z0)==n-1
                z0 = [0;z0];
            end
        end
        y0 = [log(diff(z0(1:nb)));real(z0(nb+1));log(-diff(z0(nb+1:n)))];
        
        % Find prevertices (solve param problem)
        
        % Set up normalized lengths for nonlinear equations:
        % indices of left and right integration endpoints
        left = [1,1:n-1]';
        right = [n,2:n]';
        % delete indices corresponding to vertices at Inf
        left([find(atinf)+1;nb+1]) = [];
        right([find(atinf);nb+1]) = [];
        cmplx = ((right-left) == 2);
        cmplx(1) = 0;
        cmplx(2) = 1;
        % normalize lengths
        nmlen = (w(right)-w(left))/(w(n)-w(1));
        % abs value for finite ones
        nmlen(~cmplx) = abs(nmlen(~cmplx));
        % first entry is useless (=1)
        nmlen(1) = [];
        
        % Solve nonlinear system of equations:
        
        % package data
        fdat = {n,nb,beta,nmlen,left,right,logical(cmplx),qdat};
        % set options
        opt = zeros(16,1);
        opt(1) = trace;
        opt(2) = method;
        opt(6) = 100*(n-1);
        opt(8) = tol;
        opt(9) = tol/10;
        opt(12) = nqpts;
        try
            [y,termcode] = nesolve(@stripmap.stpfun,y0,opt,fdat);
        catch
            % Have to delete the "waitbar" figure if interrupted
            close(findobj(allchild(0),'flat','Tag','TMWWaitbar'));
            error(lasterr)
        end
        if termcode~=1
            warning('Nonlinear equations solver did not terminate normally.')
        end
        
        
        % Convert y values to z
        z = zeros(n,1);
        z(2:nb) = cumsum(exp(y(1:nb-1)));
        z(nb+1:n) = i+cumsum([y(nb);-exp(y(nb+1:n-1))]);
        z = [-Inf;z(1:nb);Inf;z(nb+1:n)];
        
        % Determine multiplicative constant
        mid = mean(z(2:3));
        g = stripmap.stquad(z(3),mid,3,z,beta,qdat) - ...
            stripmap.stquad(z(2),mid,2,z,beta,qdat);
        c = (w(1)-w(2))/g;
        
        % Undo renumbering
        z(renum) = z;
        qdat(:,renum) = qdat(:,1:N);
        qdat(:,N+1+renum) = qdat(:,N+1+(1:N));
    end

    function F = stpfun(y,fdat)
        %   Returns residual for solution of nonlinear equations.
        %   $Id: stpfun.m 62 1999-01-29 00:56:34Z tad $
        
        [n,nb,beta,nmlen,left,right,cmplx,qdat] = deal(fdat{:});
        
        % In this function, n refers to the number of FINITE prevertices.
        
        % Transform y (unconstr. vars) to z (actual params)
        z = zeros(n,1);
        z(2:nb) = cumsum(exp(y(1:nb-1)));
        z(nb+1:n) = i+cumsum([y(nb);-exp(y(nb+1:n-1))]);
        
        % Compute the integrals
        zleft = z(left);
        zright = z(right);
        mid = mean([zleft.' ; zright.']).';
        c2 = cmplx;
        c2(2) = 0;
        mid(c2) = mid(c2) - sign(left(c2)-nb)*i/2;
        
        % Add ends of strip to z, and modify singularity indices
        zs = [-Inf;z(1:nb);Inf;z(nb+1:n)];
        left = left + 1 + (left > nb);
        right = right + 1 + (right > nb);
        
        % Do those staying on a side
        ints = zeros(n-1,1);
        c2(1) = 1;
        id = ~c2;
        ints(id) = stripmap.stquadh(zleft(id),mid(id),left(id),zs,beta,qdat) - ...
            stripmap.stquadh(zright(id),mid(id),right(id),zs,beta,qdat);
        
        % For the rest, go to the strip middle, across, and back to the side
        z1 = real(zleft(c2)) + i/2;
        z2 = real(zright(c2)) + i/2;
        id = ~id;
        ints(id) = stripmap.stquad(zleft(id),z1,left(id),zs,beta,qdat);
        ints(id) = ints(id) + stripmap.stquadh(z1,z2,zeros(size(z1)),zs,beta,qdat);
        ints(id) = ints(id) - stripmap.stquad(zright(id),z2,right(id),zs,beta,qdat);
        
        absval = abs(ints(~cmplx)); 		% absval(1) = abs(ints(1))
        if ~absval(1)
            rat1 = 0;
            rat2 = 0;
        else
            rat1 = absval(2:length(absval))/absval(1);
            rat2 = ints(cmplx)/ints(1);
        end
        
        if any([rat1;rat2]==0) | any(isnan([rat1;rat2])) | any(isinf([rat1;rat2]))
            % Singularities were too crowded.
            warning('Severe crowding')
        end
        
        % Compute nonlinear equation residual values.
        cmplx2 = cmplx(2:length(cmplx));
        if ~isempty(rat1)
            F1 = log( rat1 ./ nmlen(~cmplx2) );
        else
            F1 = [];
        end
        if ~isempty(rat2)
            F2 = log( rat2 ./ nmlen(cmplx2) );
        else
            F2 = [];
        end
        F = [F1;real(F2);imag(F2)];
    end
    
    function [H,RE,IM] = stplot(w,beta,z,c,re,im,options)
        %STPLOT Image of cartesian grid under Schwarz-Christoffel strip map.
        %   STPLOT(W,BETA,Z,C) will adaptively plot the images under the
        %   Schwarz-Christoffel exterior map of ten evenly spaced horizontal and
        %   vertical lines in the upper half-plane. The abscissae of the
        %   vertical lines will bracket the finite extremes of real(Z).  The
        %   arguments are as in STPARAM.
        %
        %   STPLOT(W,BETA,Z,C,M,N) will plot images of M evenly spaced vertical
        %   and N evenly spaced horizontal lines.  Horizontal lines are spaced
        %   to bracket real(Z); vertical lines are evenly spaced between 0 and
        %   1.
        %
        %   STPLOT(W,BETA,Z,C,RE,IM) will plot images of vertical lines whose
        %   real parts are given in RE and horizontal lines whose imaginary
        %   parts are given in IM.  Either argument may be empty.
        %
        %   STPLOT(W,BETA,Z,C,RE,IM,OPTIONS) allows customization of HPPLOT's
        %   behavior.  See SCPLTOPT.
        %
        %   H = STPLOT(W,BETA,Z,C,...) returns a vector of handles to all the
        %   curves drawn in the interior of the polygon.  [H,RE,IM] =
        %   STPLOT(W,BETA,Z,C,...) also returns the abscissae and ordinates of
        %   the lines comprising the grid.
        %
        %   See also SCPLTOPT, STPARAM, STMAP, STDISP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: stplot.m 298 2009-09-15 14:36:37Z driscoll $
        
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
        minre = min(real(z(~isinf(z))));
        maxre = max(real(z(~isinf(z))));
        if (length(re)==1) && (re == round(re))
            if re < 1
                re = [];
            elseif re < 2
                re = mean([minre,maxre]);
            else
                m = re;
                re = linspace(minre,maxre,m);
                dre = diff(re(1:2));
                re = linspace(minre-dre,maxre+dre,m);
            end
        end
        if (length(im)==1) && (im == round(im))
            if im < 1
                im = [];
            else
                m = im;
                im = linspace(0,1,m+2);
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
            zp = re(j) + 1i*linspace(0,1,15).';
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
                neww = stripmap.stmap(zp(new),w,beta,z,c,qdat);
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
                        'xdata',re(j)*[1 1],'ydata',[0 1])
                end
            else
                clearpoints(linh(j,1))
                addpoints(linh(j,1),real(wp(~isnan(wp))),imag(wp(~isnan(wp))));
                set(linh(j,1),'marker','none','linestyle','-','user',zp)
                if draw2
                    % Replace the points with smooth curve
                    clearpoints(linh(j,2))
                    addpoints(linh(j,2),re(j)*[1 1],[0 1])
                    set(linh(j,2),'marker','none','linestyle','-')
                end
                
            end
            drawnow
        end
        
        % Plot horizontal lines...
        x1 = min(-5,minre);
        x2 = max(5,maxre);
        wl = w(isinf(z) & z < 0);		% image of left end
        wr = w(isinf(z) & z > 0);		% image of right end
        linh1 = linh;
        linh = gobjects(length(im),2);
        for j = 1:length(im)
            % Start evenly spaced
            zp = [-Inf linspace(x1,x2,15) Inf].' + 1i*im(j);
            new = true(size(zp));
            new([1 end]) = false;
            wp = NaN(length(zp),1);
            wp([1 end]) = [wl; wr];
            
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
                neww = stripmap.stmap(zp(new),w,beta,z,c,qdat);
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
                        'xdata',real(zp([2 end-1])),'ydata',im(j)*[1 1])
                end
            else
                clearpoints(linh(j,1))
                addpoints(linh(j,1),real(wp(~isnan(wp))),imag(wp(~isnan(wp))));
                set(linh(j,1),'marker','none','linestyle','-','user',zp)
                if draw2
                    % Replace the points with (hopefully) a smooth circle
                    clearpoints(linh(j,2))
                    addpoints(linh(j,2),real(zp([2 end-1])),im(j)*[1 1])
                    set(linh(j,2),'marker','none','linestyle','-')
                end
                
            end
            drawnow
        end
        
        % Finish up
        
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
    
    function I = stquad(z1,z2,sing1,z,beta,qdat)
        %STQUAD (not intended for calling directly by the user)
        %   Numerical quadrature for the strip map.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: stquad.m 298 2009-09-15 14:36:37Z driscoll $
        
        %   z1,z2 are vectors of left and right endpoints.  sing1 is a vector of
        %   integer indices which label the singularities in z1.  So if sing1(5)
        %   = 3, then z1(5) = z(3).  A zero means no singularity.  z is the
        %   vector of *all* singularities, including the "ends" of the strip at
        %   \pm Inf. beta is the vector of associated turning angles. qdat is
        %   quadrature data from SCQDATA. It should include all the beta values,
        %   even though the ends are never used in this manner.
        %
        %   Make sure z and beta are column vectors.
        %
        %   STQUAD integrates from a possible singularity at the left end to a
        %   regular point at the right.  If both endpoints are singularities,
        %   you must break the integral into two pieces and make two calls.
        %
        %   The integral is subdivided, if necessary, so that no singularity
        %   lies closer to the left endpoint than 1/2 the length of the
        %   integration (sub)interval.
        
        n = length(z);
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
            ind = rem(sng+n,n+1)+1;
            % Adjust Gauss-Jacobi nodes and weights to interval.
            nd = ((zr-za)*qdat(:,ind) + zr + za).'/2; % G-J nodes
            wt = ((zr-za)/2) * qdat(:,ind+n+1); 	% G-J weights
            if any( diff([za;nd(:);zr])==0 )
                % Endpoints are practically coincident.
                I(k) = 0;
            else
                % Use Gauss-Jacobi on first subinterval, if necessary.
                if sng > 0
                    wt = wt*(abs(zr-za)/2)^beta(sng);
                end
                I(k) = stripmap.stderiv(nd,z,beta,1,sng)*wt;
                while (dist < 1) & ~isnan(I(k))
                    % Do regular Gaussian quad on other subintervals.
                    zl = zr;
                    dist = min(1,2*min(abs(z-zl))/abs(zl-zb));
                    zr = zl + dist*(zb-zl);
                    nd = ((zr-zl)*qdat(:,n+1) + zr + zl).'/2;
                    wt = ((zr-zl)/2) * qdat(:,2*n+2);
                    I(k) = I(k) + stripmap.stderiv(nd,z,beta,1)*wt;
                end
            end
        end
    end
    
    function I = stquadh(z1,z2,sing1,z,beta,qdat)
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: stquadh.m 298 2009-09-15 14:36:37Z driscoll $
        
        %   STQUAD applies the "1/2 rule" by assuming that the distance from the
        %   integration interval to the nearest singularity is equal to the
        %   distance from the left endpoint to the nearest singularity. This is
        %   certainly true e.g. if one begins integration at the nearest
        %   singularity to the target point. However, it may be violated in
        %   important circumstances, such as when one integrates from the
        %   next-nearest singularity (because the nearest maps to inf), or when
        %   one is integrating between adjacent prevertices (as in the param
        %   problem). The main difficulty is the possibility of singularities
        %   from the "other side" of the strip intruding.
        %
        %   Here we assume that the integration intervals are horizontal. This
        %   function recursively subdivides the interval until the "1/2 rule" is
        %   satisfied for singularities "between" the endpoints. Actually, we
        %   use a more stringent "alpha rule", for alpha > 1/2, as this seems to
        %   be necessary sometimes.
        %
        %   There must be no singularities *inside* the interval, of course.
        
        n = length(z);
        if isempty(sing1)
            sing1 = zeros(length(z1),1);
        end
        
        I = zeros(size(z1));
        nontriv = find(z1(:)~=z2(:))';
        
        for k = nontriv
            za = z1(k);
            zb = z2(k);
            sng = sing1(k);
            
            % alf==1/2 means the "1/2 rule." Better to be more strict.
            alf = .75;
            
            % Given integration length
            d = real(zb) - real(za);
            
            % Compute horizontal position (signed) and vertical distance (positive)
            % from the singularities to the left endpoint. If we are going from right
            % to left, reverse the sense of horizontal.
            dx = (real(z) - real(za)) * sign(d);
            dy = abs(imag(z) - imag(za));
            
            % We have to be concerned with singularities lying to the right (left if
            % d<0) of the left integration endpoint.
            toright = (dx > 0) & (~isinf(z));
            % For points with small enough dx, the limitation is purely due to dy. For
            % others, it must be calculated.
            active = ( dx > dy/alf ) & toright;
            % Make sure that the left endpoint won't be included
            if sng
                active(sng) = 0;
            end
            
            % For those active, find the integration length constraint. This comes
            % from making the sing/right-endpoint distance equal to alf*L.
            x = dx(active);
            y = dy(active);
            L = (x - sqrt((alf*x).^2 - (1-alf^2)*y.^2)) / (1-alf^2);
            
            % What is the maximum allowable integration length?
            L = min([ L(:); dy(toright & ~active)/alf ]);
            
            if L < abs(d)
                % Apply STQUAD on the safe part and recurse on the rest
                zmid = za + L*sign(d);
                I(k) = stripmap.stquad(za,zmid,sng,z,beta,qdat);
                I(k) = I(k) + stripmap.stquadh(zmid,zb,0,z,beta,qdat);
            else
                % No restriction
                I(k) = stripmap.stquad(za,zb,sng,z,beta,qdat);
            end
            
        end
    end
end

end
