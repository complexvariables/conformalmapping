classdef diskmap < scmap
%DISKMAP Schwarz-Christoffel disk map object.
%   DISKMAP(P) constructs a Schwarz-Christoffel disk map object for the
%   polygon P. The parameter problem is solved using default options for
%   the prevertices and the multiplicative constant.
%
%   DISKMAP(P,OPTIONS) uses an options structure of the type created by
%   SCMAPOPT in solving the parameter problem.
%
%   DISKMAP(P,Z) creates a diskmap object having the given prevertices Z
%   (the mulitiplicative constant is found automatically). There is no
%   checking to ensure that the prevertices are consistent with the
%   given polygon. DISKMAP(P,Z,C) also uses the supplied constant. An
%   OPTIONS argument can be added, although only the error tolerance
%   will be used.
%
%   DISKMAP(M), where M is a diskmap object, just returns M.
%
%   DISKMAP(M,P) returns a new diskmap object for the polygon P using
%   the options in diskmap M. The prevertices of M will be used as the
%   starting guess for the parameter problem of the new map. Thus P
%   should properly be a perturbation of the polygon for M. An OPTIONS
%   structure may be added to override options in M.
%
%   DISKMAP(Z,ALPHA) creates a map using the given prevertices and the
%   interior polygon angles described by ALPHA (see POLYGON help). The
%   image polygon is deduced by computing S-C integrals assuming a
%   multiplicative constant of 1. DISKMAP(Z,ALPHA,C) uses the given
%   constant instead.
%
%   See also SCMAPOPT, classes POLYGON, SCMAP.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll.

properties(SetAccess=protected)
    prevertex
    constant
    qdata
    accuracyVal
    centerVal
end

methods
    function map = diskmap(varargin)
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
            case 'diskmap'
                oldmap = varargin{1};
                % Continuation of given map to given polygon
                poly = varargin{2};
                opt = scmapopt(oldmap);
                z0 = prevertex(oldmap);
                if length(z0) ~= length(poly)
                    msg = 'Polygon %s must have the same length as that in %s.';
                    error(sprintf(msg,inputname(2),inputname(1)))
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
                        error(sprintf(msg,inputname(j)))
                    end
                end
                
            case 'double'
                % Args are the prevertex vector, then angle vector
                z = varargin{1};
                alpha = varargin{2};
                poly = polygon(NaN*alpha*1i,alpha);
                c = 1;
                for j = 3:length(varargin)
                    if isa(varargin{j},'struct')
                        opt = varargin{j};
                    elseif length(varargin{j})==1
                        c = varargin{j};
                    else
                        msg = 'Unable to parse argument ''%s''.';
                        error(sprintf(msg,inputname(j)))
                    end
                end
                
            otherwise
                msg = 'Expected ''%s'' to be a polygon, a diskmap, or prevertices.';
                error(sprintf(msg,inputname(1)))
        end % input switch
        
        % Retrieve options
        opt = scmapopt(opt);
        
        % Take actions based on what needs to be filled in
        
        if isempty(z)
            % Solve parameter problem
            % Apply SCFIX to enforce solver rules
            [w,beta] = scfix('d',vertex(poly),angle(poly)-1);
            poly = polygon(w,beta+1);             % in case polygon was altered
            
            [z,c,qdata] = diskmap.dparam(w,beta,opt.InitialGuess,opt);
        end
        
        if isempty(qdata)
            % Base accuracy of quadrature on given options
            nqpts = ceil(-log10(opt.Tolerance));
            qdata = scqdata(angle(poly)-1,nqpts);
        end
        
        if isempty(c)
            % Find constant
            w = vertex(poly);
            beta = angle(poly)-1;
            idx = 1 + find(~isinf(z(2:end)), 1 );
            mid = mean(z([1 idx]));
            I = diskmap.dquad(z(1),mid,1,z,beta,qdata) ...
                - diskmap.dquad(z(idx),mid,idx,z,beta,qdata);
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
        
        % Find conformal center
        map.centerVal = center(map);
        
        % Fill in apparent accuracy
        map.accuracyVal = accuracy(map);        
    end
    
    function acc = accuracy(M)
        %ACCURACY Apparent accuracy of Schwarz-Christoffel disk map.
        %   ACCURACY(M) estimates the accuracy of the Schwarz-Christoffel disk
        %   map M. The technique used is to compare the differences between
        %   successive finite vertices to the integral between the corresponding
        %   prevertices, and return the maximum.
        %
        %   See also DISKMAP.
        
        % If an accuracy has been assigned, don't question it
        if ~isempty(M.accuracyVal)
            acc = M.accuracyVal;
            return
        end
        
        % Get data for low-level functions
        p = M.polygon;
        w = vertex(p);
        beta = angle(p) - 1;
        z = M.prevertex;
        c = M.constant;
        qdata = M.qdata;
        
        % Test accuracy by integrating between consecutive finite prevertices, and
        % comparing to differences of vertices.
        
        n = length(w);
        idx = find(~isinf(w));
        wf = w(idx);				% finite vertices
        
        % Two columns hold endpoint indices for integrations
        idx = [idx(1:end) idx([2:end 1])];
        
        % Always use center as the integration midpoint
        %dtheta = mod(angle(z(idx(:,2))./z(idx(:,1))),2*pi);
        %mid = z(idx(:,1)).*exp(i*dtheta/2);
        mid = zeros(length(idx),1);
        
        % Do the integrations
        I = M.dquad(z(idx(:,1)),mid,idx(:,1),z,beta,qdata) - ...
            M.dquad(z(idx(:,2)),mid,idx(:,2),z,beta,qdata);
        
        acc = max(abs( c*I - diff(wf([1:end 1])) ));
    end
    
    function wc = center(map,wc)
        %CENTER Conformal center of Schwarz-Christoffel disk map.
        %   CENTER(M) returns the conformal center (image of 0) of the
        %   Schwarz-Christoffel disk map represented by M.
        %
        %   CENTER(M,WC) computes a map conformally equivalent to M but with
        %   conformal center WC (provided WC is inside the polygon of M), and
        %   returns the new map. If WC is empty, you will be asked to select it
        %   graphically.
        %
        %   See also DISKMAP.
        
        if nargin == 1
            % Return center
            wc = map.centerVal;
            if isempty(wc)
                p = map.polygon;
                wc = map.dmap(0,vertex(p),angle(p)-1,...
                    map.prevertex,map.constant,map.qdata);
            end
            
        else
            % Set center
            p = map.polygon;
            qdata = map.qdata;
            z = map.prevertex;
            w = vertex(p);
            beta = angle(p) - 1;
            
            if isempty(wc)
                fig = figure;
                plot(p)
                title('Click at conformal center')
                [xc,yc] = ginput(1);
                wc = xc + 1i*yc;
                delete(fig)
            end
            
            if ~any(isinf(w)) && ~isinpoly(wc,p)
                error('Conformal center must be inside polygon.')
            end
            
            % Find inverse image of wc under current map
            zc = map.dinvmap(wc,w,beta,z,map.constant,qdata);
            
            % Use Moebius transform to reset prevertices
            y = ((1-zc')/(1-zc))*(z-zc)./(1-zc'*z);
            y(length(y)) = 1;			% force it to be exact
            y = y./abs(y);
            
            % Recalculate constant
            mid = mean(y(1:2));
            I = map.dquad(y(1),mid,1,y,beta,qdata) ...
                - map.dquad(y(2),mid,2,y,beta,qdata);
            c = diff(w(1:2))/I;
            
            % Assign new values
            map.prevertex = y;
            map.constant = c;
            map.centerVal = wc;
            map.accuracyVal = [];
            map.accuracyVal = accuracy(map);
            wc = map;
        end
    end
    
    function out = char(f)
        %CHAR   Pretty-print a Schwarz-Christoffel disk map.
        
        p = f.polygon;
        w = vertex(p);
        alpha = angle(p);
        z = f.prevertex;
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
        L{end+1}=sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c)));
        
        wc = center(f);
        if imag(wc) < 0
            s = '-';
        else
            s = '+';
        end
        L{end+1}=sprintf('  Conformal center at %.4f %c %.4fi',real(wc),s,abs(imag(wc)));
        
        L{end+1} = sprintf('  Apparent accuracy is %.2e',f.accuracyVal);
        L{end+1} = ' ';
        
        out = L;
    end
    
    function out = disp(M)
        %DISPLAY Display parameters of a Schwarz-Christoffel disk map.
        
        p = M.polygon;
        w = vertex(p);
        alpha = angle(p);
        z = M.prevertex;
        c = M.constant;
        
        L = {' '; '  diskmap object:';  ' '};
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
        L{end+1}=sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c)));
        
        wc = center(M);
        if imag(wc) < 0
            s = '-';
        else
            s = '+';
        end
        L{end+1}=sprintf('  Conformal center at %.4f %c %.4fi',real(wc),s,abs(imag(wc)));
        
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
        %EVAL Evaluate Schwarz-Christoffel disk map at points.
        %   EVAL(M,ZP) evaluates the Schwarz-Christoffel map M at the points ZP
        %   in the unit disk. The default tolerance of M is used.
        %
        %   EVAL(M,ZP,TOL) attempts to give an answer accurate to TOL. If TOL is
        %   less than the accuracy of M, this is unlikely to be met.
        %
        %   See also DISKMAP, EVALINV.
        
        if nargin < 3
            qdata = M.qdata;
        else
            qdata = tol;
        end
        
        p = M.polygon;
        wp = NaN*zp;
        idx = abs(zp) <= 1+eps;
        wp(idx) = ...
            M.dmap(zp(idx),vertex(p),angle(p)-1,M.prevertex,M.constant,qdata);
    end
    
    function fp = evaldiff(M,zp)
        %EVALDIFF Derivative of Schwarz-Christoffel disk map at points.
        %   EVALDIFF(M,ZP) computes the derivative of the Schwarz-Christoffel
        %   disk map M at the points ZP.
        %
        %   See also DISKMAP, EVAL.
        
        z = M.prevertex;
        c = M.constant;
        beta = angle(M.polygon) - 1;
        
        fp = M.dderiv(zp,z,beta,c);
    end
    
    function [zp,flag] = evalinv(M,wp,tol,z0)
        %EVALINV Invert Schwarz-Christoffel disk map at points.
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
        %   [ZP,FLAG] = EVALINV(...) also returns a vector of indices where the
        %   method was unable to produce a sufficiently small residual. A warning
        %   is issued when this occurs.
        %
        %   See also DISKMAP, EVAL.
        
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
        [zp,flag] = dinvmap(wp,vertex(p),angle(p)-1,M.prevertex,M.constant,...
            qdata,z0,[0 tol]);
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
        %   Given a diskmap M, FORWARDPOLY(M) returns the polygon that is
        %   formed using the prevertices, angles, and quadrature data of that
        %   map. If the prevertices were found from the solution of a
        %   parameter problem, then the result should agree closely with the
        %   original polygon that was supplied.
        
        z = map.prevertex;
        alpha = angle(map.polygon);
        c = map.constant;
        
        n = length(z);
        
        % Since there is no parameter problem, use high accuracy in quadrature.
        qdata = scqdata(alpha-1,16);
        
        w = zeros(n,1);
        atinf = (alpha < eps);
        w(atinf) = Inf;
        
        % Endpoints of integrations
        idx = find(~atinf);
        idx = [idx(1:end-1) idx(2:end)];
        
        % Origin is midpoint of every integration
        mid = zeros(length(idx),1);
        
        % Integrations
        I = map.dquad(z(idx(:,1)),mid,idx(:,1),z,alpha-1,qdata) - ...
            map.dquad(z(idx(:,2)),mid,idx(:,2),z,alpha-1,qdata);
        
        % Deduce vertices
        w(~atinf) = c*cumsum([0;I]);
        
        p = polygon(w,alpha);
    end
    
    function varargout = get(map,varargin)
        %GET    Get map parameters.
        %   [VAL1,VAL2,...] = GET(F,'PROP1','PROP2',...) returns the values of the
        %   map F corresponding to the requested properties. Valid properties
        %   are:
        %
        %       polygon, options, prevertex, constant, center
        
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
                case 'cen'
                    varargout{j} = center(map);
                otherwise
                    warning(sprintf('Property ''%s'' not recognized.\n',varargin{j}))
                    varargout{j} = [];
            end
        end
    end
    
    function M1 = hplmap(M)
        %HPLMAP Convert Schwarz-Christoffel disk map to a map from the half-plane.
        
        p = M.polygon;
        [z1,c1] = M.disk2hp(vertex(p),angle(p)-1,M.prevertex,M.constant);
        M1 = hplmap(p,scmapopt(M),z1,c1);
    end
    
    function M = mtimes(M,c)
        %   Scale the image of a map by a complex constant.
        
        % May need to swap arguments
        if isa(M,'double') & isa(c,'diskmap')
            tmp = M;
            M = c;
            c = tmp;
        end
        
        M.constant = c*M.constant;
        M.polygon = c*M.polygon;
        M.centerVal = c*M.centerVal;
    end
    
    function v = parameters(M)
        %PARAMETERS Return a structure of the Schwarz-Christoffel map parameters.
        
        v.prevertex = M.prevertex;
        v.constant = M.constant;
    end
    
    function [h,r,theta] = plot(M,varargin)
        %PLOT Visualize a Schwarz-Christoffel disk map.
        %   PLOT(M) plots the polygon associated with the Schwarz-Christoffel
        %   disk map M and the images of ten evenly spaced circles and radii
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
        %   See also DISKMAP, EVAL.
        
        p = M.polygon;
        w = vertex(p);
        beta = angle(p) - 1;
        z = M.prevertex;
        c = M.constant;
        
        if nargin == 1
            [a1,a2,a3] = M.dplot(w,beta,z,c);
        elseif length(varargin) == 1
            % Tolerance given only
            [a1,a2,a3] = M.dplot(w,beta,z,c,10,10,ceil(-log10(varargin{1})));
        elseif length(varargin) == 2
            % R, theta given only
            [a1,a2,a3] = M.dplot(w,beta,z,c,varargin{1},varargin{2});
        else
            % All given
            nqpts = ceil(-log10(varargin{3}));
            [a1,a2,a3] = M.dplot(w,beta,z,c,varargin{1},varargin{2},nqpts);
        end
        
        if nargout > 0
            h = a1;
            r = a2;
            theta = a3;
        end
    end
    
    function M = plus(M,a)
        %   Add a constant to the image of a diskmap (i.e., translate image).
        
        % May need to swap arguments
        if isa(M,'double') & isa(a,'diskmap')
            tmp = M;
            M = a;
            a = tmp;
        end
        
        if length(a)==1 & isa(a,'double')
            M.centerVal = M.centerVal + a;
            M.polygon = M.polygon + a;
        else
            error('Addition is not defined for these operands.')
        end
    end
end

methods(Hidden,Static)
    function fprime = dderiv(zp,z,beta,c)
        %DDERIV Derivative of the disk map.
        %   DDERIV(ZP,Z,BETA,C) returns the derivative at the points of ZP of
        %   the Schwarz-Christoffel disk map defined by Z, BETA, and C.
        %
        %   See also DPARAM, DMAP.
        
        % Support old syntax
        if nargin < 4
            c = 1;
        end
        
        z = z(:);
        beta = beta(:);
        zprow = zp(:).';
        fprime = zeros(size(zp));
        
        npts = length(zp(:));
        terms = 1 - zprow(ones(length(beta),1),:)./z(:,ones(npts,1));
        fprime(:) = c*exp(sum(log(terms).*beta(:,ones(npts,1))));
    end
    
    function ddisp(w,beta,z,c)
        %DDISP  Display results of Schwarz-Christoffel disk parameter problem.
        %   DDISP(W,BETA,Z,C) displays the results of DPARAM in a pleasant way.
        %
        %   See also DPARAM, DPLOT.
        
        disp(' ')
        disp('      vertex [w]          beta        prevertex [z]         arg(z)/pi')
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
        disp(sprintf('  c = %.8g %c %.8gi\n',real(c),s,abs(imag(c))))
    end
    
    function [y,d] = dfixwc(w,beta,z,c,wc,tol)
        %DFIXWC Fix conformal center of disk map.
        %   The conformal center WC of a Schwarz-Christoffel interior disk map
        %   is defined as the image of zero.  The parameter problem solver
        %   DPARAM does not allow control over the placement of the conformal
        %   center.  Using the output Z,C from DPARAM, [Z0,C0] =
        %   DFIXWC(W,BETA,Z,C,WC) computes a Moebius transformation so that if
        %   Z0 and C0 are used in place of Z and C, the conformal center of the
        %   resulting map will be WC.
        %
        %   [Z0,C0] = DFIXWC(W,BETA,Z,C,WC,TOL) uses tolerance TOL.
        %
        %   See also DPARAM, PTSOURCE.
        
        n = length(w);
        
        if nargin < 6
            [trace,tol,method] = sctool.scparopt([]);
        end
        
        zc = dinvmap(wc,w,beta,z,c,tol);
        
        % Transform prevertices.
        y = ((1-zc')/(1-zc))*(z-zc)./(1-zc'*z);
        y(n) = 1;				% force it to be exact
        y = y./abs(y);
        
        % Recalculate constant from scratch.
        mid = (y(1)+y(2))/2;
        qdat = scqdata(beta,ceil(-log10(tol)));
        d = (w(1) - w(2))/...
            (diskmap.dquad(y(2),mid,2,y,beta,qdat) ...
            - diskmap.dquad(y(1),mid,1,y,beta,qdat));
    end
    
    function zdot = dimapfun(wp,yp,scale,z,beta,c)
        %   Used by DINVMAP for solution of an ODE.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: dimapfun.m 7 1998-05-10 04:37:19Z tad $
        
        lenyp = length(yp);
        lenzp = lenyp/2;
        zp = (yp(1:lenzp)+sqrt(-1)*yp(lenzp+1:lenyp));
        
        f = scale./diskmap.dderiv(zp,z,beta,c);
        zdot = [real(f);imag(f)];        
    end
    
    function [zp,flag] = dinvmap(wp,w,beta,z,c,qdat,z0,options)
        %DINVMAP Schwarz-Christoffel disk inverse map.
        %   DINVMAP(WP,W,BETA,Z,C,TOL) computes the inverse of the
        %   Schwarz-Christoffel disk map (i.e., from a polygon to the disk) at
        %   the points given in vector WP. The other arguments are as in
        %   DPARAM. TOL is a scalar tolerance, or a quadrature-data matrix QDAT
        %   as returned by SCQDATA, or may be omitted.
        %
        %   The default algorithm is to solve an ODE in order to obtain a fair
        %   approximation for ZP, and then improve ZP with Newton iterations.
        %   The ODE solution at WP requires a vector Z0 whose forward image W0
        %   is such that for each j, the line segment connecting WP(j) and W0(j)
        %   lies inside the polygon.  By default Z0 is chosen by a fairly robust
        %   automatic process.  Using a parameter (see below), you can choose to
        %   use either an ODE solution or Newton iterations exclusively.
        %
        %   DINVMAP(WP,W,BETA,Z,C,TOL,Z0) has two interpretations.  If the ODE
        %   solution is being used, Z0 overrides the automatic selection of
        %   initial points.  (This can be handy in convex polygons, where the
        %   choice of Z0 is trivial.)  Otherwise, Z0 is taken as an initial
        %   guess to ZP.  In either case, if length(Z0)==1, the value Z0 is used
        %   for all elements of WP; otherwise, length(Z0) should equal
        %   length(WP).
        %
        %   DINVMAP(WP,W,BETA,Z,C,TOL,Z0,OPTIONS) uses a vector of parameters
        %   that control the algorithm.  See SCINVOPT.
        %
        %   [ZP,FLAG] = DINVMAP(...) also returns a vector of indices where the
        %   method was unable to produce a sufficiently small residual. A warning
        %   is issued when this occurs.
        %
        %   See also SCINVOPT, DPARAM, DMAP.
        
        n = length(w);
        beta = beta(:);
        z = z(:);
        zp = zeros(size(wp));
        wp = wp(:);
        lenwp = length(wp);
        
        if nargin < 8
            options = [];
            if nargin < 7
                z0 = [];
                if nargin < 6
                    qdat = [];
                end
            end
        end
        
        [ode,newton,tol,maxiter] = sctool.scinvopt(options);
        
        if isempty(qdat)
            qdat = tol;
        end
        
        if length(qdat)==1
            qdat = sctool.scqdata(beta,max(ceil(-log10(qdat)),2));
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
        if lenwp==0, flag = []; return, end
        
        % ODE
        if ode
            if isempty(z0)
                % Pick a value z0 (not a singularity) and compute the map there.
                map = @(zp) diskmap.dmap(zp,w,beta,z,c,qdat);
                [z0,w0] = sctool.findz0('d',wp(~done),map,w,beta,z,c,qdat);
            else
                w0 = diskmap.dmap(z0,w,beta,z,c,qdat);
                if length(z0)==1 && lenwp > 1
                    z0 = z0(:,ones(lenwp,1)).';
                    w0 = w0(:,ones(lenwp,1)).';
                end
                w0 = w0(~done);
                z0 = z0(~done);
            end
            
            % Use relaxed ODE tol if improving with Newton.
            odetol = max(tol,1e-4*(newton));
            opt = odeset('abstol',odetol,'reltol',odetol);
            
            % Rescale dependent coordinate
            scale = (wp(~done) - w0(:));
            
            % Solve ODE
            z0 = [real(z0);imag(z0)];
            odefun = @(w,y) diskmap.dimapfun(w,y,scale,z,beta,c);
            [t,y] = ode113(odefun,[0,0.5,1],z0,opt);
            [m,leny] = size(y);
            zp(~done) = y(m,1:lenwp)+sqrt(-1)*y(m,lenwp+1:leny);
            abszp = abs(zp);
            out = abszp > 1;
            zp(out) = zp(out)./abszp(out);
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
            while ~all(done) && k < maxiter
                F = wp(~done) - diskmap.dmap(zn(~done),w,beta,z,c,qdat);
                m = length(F);
                dF = c*exp(sum(beta(:,ones(m,1)).*...
                    log(1-(zn(~done,ones(n,1)).')./z(:,ones(m,1)))));
                zn(~done) = zn(~done) + F(:)./dF(:);
                out = abs(zn) > 1;
                zn(out) = sign(zn(out));
                done(~done) = (abs(F)< tol);
                k = k+1;
            end
            if any(abs(F)> tol)
                str = sprintf('Check solution; maximum residual = %.3g\n',max(abs(F)));
                warning(str)
            end
            zp(:) = zn;
        end;
        
        flag = find(~done);
        
    end
    
    function [zhp,chp] = disk2hp(w,beta,z,c)
        %DISK2HP Convert solution from the disk to one from the half-plane.
        %   [ZHP,CHP] = DISK2HP(W,BETA,Z,C) quickly transforms the solution Z,C
        %   of the Schwarz-Christoffel disk mapping parameter problem to the
        %   solution ZHP,CHP of the half-plane problem.
        %
        %   See also HP2DISK, DPARAM, HPPARAM.
        
        n = length(w);
        zhp = zeros(size(z));
        zhp(n) = Inf;
        zhp(1:n-1) = -i*(z(1:n-1)+1)./(z(1:n-1)-1); % Mobius transfmn
        zhp = real(zhp);
        
        % Recalculate constant from scratch.
        mid = mean(zhp(1:2));
        qdat = sctool.scqdata(beta(1:n-1),16);
        chp = (w(1)-w(2))/(hplmap.hpquad(zhp(2),mid,2,zhp(1:n-1),beta(1:n-1),qdat) - ...
            hplmap.hpquad(zhp(1),mid,1,zhp(1:n-1),beta(1:n-1),qdat));
    end
    
    function wp = dmap(zp,w,beta,z,c,qdat)
        %DMAP  Schwarz-Christoffel disk map.
        %   DMAP(ZP,W,BETA,Z,C,QDAT) computes the values of the Schwarz-
        %   Christoffel disk map at the points in vector ZP. The arguments W,
        %   BETA, Z, C, and QDAT are as in DPARAM. DMAP returns a vector the
        %   same size as ZP.
        %
        %   DMAP(ZP,W,BETA,Z,C,TOL) uses quadrature data intended to give an
        %   answer accurate to within TOL.
        %
        %   DMAP(ZP,W,BETA,Z,C) uses a tolerance of 1e-8.
        %
        %   See also DPARAM, DPLOT, DINVMAP.
        
        if isempty(zp)
            wp = [];
            return
        end
        
        n = length(z);
        w = w(:);
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
        
        % "Bad" points are closest to a prevertex of infinity.
        atinf = find(isinf(w)); 		% infinite vertices
        bad = ismember(sing,atinf) & ~vertex;
        
        if any(bad)
            % Can't integrate starting at pre-infinity: find conformal center to use
            % as integration basis.
            if ~isinf(w(n-1))
                wc = w(n-1) + c*diskmap.dquad(z(n-1),0,n-1,z,beta,qdat);
            else
                wc = w(n) + c*diskmap.dquad(z(n),0,n,z,beta,qdat);
            end
        end
        
        % zs = the starting singularities
        zs = z(sing);
        % ws = map(zs)
        ws = w(sing);
        
        % Compute the map directly at "normal" points.
        normal = ~bad & ~vertex;
        if any(normal)
            I = diskmap.dquad(zs(normal),zp(normal),sing(normal),z,beta,qdat);
            wp(normal) = ws(normal) + c*I;
        end
        
        % Compute map at "bad" points, using conformal center as basis, to avoid
        % integration where right endpoint is too close to a singularity.
        if any(bad)
            I = diskmap.dquad(zp(bad),zeros(sum(bad),1),...
                zeros(sum(bad),1),z,beta,qdat);
            wp(bad) = wc - c*I;
        end
        
        wp = reshape(wp,shape);
    end
    
    function [z,c,qdat] = dparam(w,beta,z0,options)
        %DPARAM Schwarz-Christoffel disk parameter problem.
        %   [Z,C,QDAT] = DPARAM(W,BETA) solves the Schwarz-Christoffel mapping
        %   parameter problem with the disk as fundamental domain and the
        %   polygon specified by W as the target. W must be a vector of the
        %   vertices of the polygon, specified in counterclockwise order, and
        %   BETA should be a vector of the turning angles of the polygon; see
        %   SCANGLE for details. If successful, DPARAM will return Z, a vector
        %   of the pre-images of W; C, the multiplicative constant of the
        %   conformal map; and QDAT, an optional matrix of quadrature data used
        %   by some of the other S-C routines.
        %
        %   [Z,C,QDAT] = DPARAM(W,BETA,Z0) uses Z0 as an initial guess for Z.
        %
        %   [Z,C,QDAT] = DPARAM(W,BETA,TOL) attempts to find an answer within
        %   tolerance TOL. (Also see SCPAROPT.)
        %
        %   [Z,C,QDAT] = DPARAM(W,BETA,Z0,OPTIONS) uses a vector of control
        %   parameters. See SCPAROPT.
        %
        %   See also SCPAROPT, DRAWPOLY, DDISP, DPLOT, DMAP, DINVMAP.
        
        %   Copyright 1998-2001 by Toby Driscoll.
        %   $Id: dparam.m 199 2002-09-13 18:54:27Z driscoll $
        
        import sctool.*
        
        n = length(w);				% no. of vertices
        w = w(:);
        beta = beta(:);
        
        % Set up defaults for missing args
        if nargin < 4
            options = [];
            if nargin < 3
                z0 = [];
            end
        end
        
        err = sccheck('d',w,beta);
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
        qdat = scqdata(beta,nqpts); 		% quadrature data
        
        atinf = (beta <= -1);
        
        if n==3
            % Trivial solution
            z = [-i;(1-i)/sqrt(2);1];
            
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
                y0 = zeros(n-3,1);
            else
                z0 = z0(:)./abs(z0(:));
                % Moebius to make th(n-2:n)=[1,1.5,2]*pi;
                Am = moebius(z0(n-2:n),[-1;-i;1]);
                z0 = Am(z0);
                th = angle(z0);
                th(th<=0) = th(th<=0) + 2*pi;
                dt = diff([0;th(1:n-2)]);
                y0 = log(dt(1:n-3)./dt(2:n-2));
            end
            
            % Solve nonlinear system of equations:
            
            % package data
            fdat = {n,beta,nmlen,left,right,logical(cmplx),qdat};
            % set options
            opt = zeros(16,1);
            opt(1) = trace;
            opt(2) = method;
            opt(6) = 100*(n-3);
            opt(8) = tol;
            opt(9) = min(eps^(2/3),tol/10);
            opt(12) = nqpts;
            try
                [y,termcode] = nesolve(@diskmap.dpfun,y0,opt,fdat);
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
            theta = pi*cs(1:n-3)/cs(n-2);
            z = ones(n,1);
            z([1:n-3]) = exp(i*theta);
            z(n-2:n-1) = [-1;-i];
        end
        
        % Determine scaling constant
        mid = (z(1)+z(2))/2;
        c = (w(1) - w(2))/...
            (diskmap.dquad(z(2),mid,2,z,beta,qdat) ...
            - diskmap.dquad(z(1),mid,1,z,beta,qdat));
        
    end
    
    function F = dpfun(y,fdat)
        %   Returns residual for solution of nonlinear equations.
        
        [n,beta,nmlen,left,right,cmplx,qdat] = deal(fdat{:});
        
        % Convert y values to z (prevertices)
        cs = cumsum(cumprod([1;exp(-y)]));
        theta = pi*cs(1:n-3)/cs(length(cs));
        z = ones(n,1);
        z(1:n-3) = exp(i*theta);
        z(n-2:n-1) = [-1;-i];
        
%         %%% Check crowding.
%         %%if any(diff(theta)<eps) | any(isnan(theta))
%         %%  % Since abs(y) is large, use it as the penalty function.
%         %%  F = y;
%         %%  disp('Warning: Severe crowding')
%         %%  return
%         %%end
        
        % Compute the integrals
        zleft = z(left);
        zright = z(right);
        angl = angle(zleft);
        mid = exp(i*(angl + rem(angle(zright./zleft)+2*pi,2*pi)/2));
        % For integrals between nonadjacent singularities, choose 0 as intermediate
        % integration point.
        mid(cmplx) = zeros(size(mid(cmplx)));
        % If any are complex, the first one must be too.
        if any(cmplx)
            cmplx(1) = 1;
        end
        
        ints = NaN*zleft;
        ints(~cmplx) = sctool.dabsquad(zleft(~cmplx),mid(~cmplx),left(~cmplx),...
            z,beta,qdat) + ...
            sctool.dabsquad(zright(~cmplx),mid(~cmplx),right(~cmplx),z,beta,qdat);
        if any(cmplx)
            ints(cmplx) = sctool.dquad(zleft(cmplx),mid(cmplx),left(cmplx),z,beta,qdat) - ...
                sctool.dquad(zright(cmplx),mid(cmplx),right(cmplx),z,beta,qdat);
        end
        
        if any(ints==0)
            % Singularities were too crowded in practice.
%             %%  F = y;
            warning('Severe crowding')
        end
        
        % Compute nonlinear equation residual values.
        cmplx(1) = 0;
        F1 = ints(~cmplx);
        F1 = F1(2:end)/abs(F1(1));
        F2 = ints(cmplx)/ints(1);
        F = [F1;real(F2);imag(F2)] - nmlen;
    end
    
    function [H,R2,THETA] = dplot(w,beta,z,c,R,theta,options)
        %DPLOT  Image of polar grid under Schwarz-Christoffel disk map.
        %   DPLOT(W,BETA,Z,C) will adaptively plot the images under the
        %   Schwarz-Christoffel disk map of ten evenly spaced circles and rays
        %   in the unit disk.  The arguments are as in DPARAM.
        %
        %   DPLOT(W,BETA,Z,C,M,N) will plot images of M evenly spaced circles
        %   and N evenly spaced rays.
        %
        %   DPLOT(W,BETA,Z,C,R,THETA) will plot images of circles whose radii
        %   are given in R and rays whose arguments are given in THETA.  Either
        %   argument may be empty.
        %
        %   DPLOT(W,BETA,Z,C,R,THETA,OPTIONS) allows customization of DPLOT's
        %   behavior.  See SCPLTOPT.
        %
        %   H = DPLOT(W,BETA,Z,C,...) returns a vector of handles to all the
        %   curves drawn in the interior of the polygon.  [H,R,THETA] =
        %   DPLOT(W,BETA,Z,C,...) also returns the moduli and arguments of the
        %   curves comprising the grid.
        %
        %   See also SCPLTOPT, DPARAM, DMAP, DDISP.
        
        %   Copyright 1998 by Toby Driscoll.
        %   $Id: dplot.m 226 2003-01-08 16:59:17Z driscoll $
        
        import sctool.*
        w = w(:);
        beta = beta(:);
        z = z(:);
        
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
            R = linspace(0,1,m);
            R([1,m]) = [];
        end
        if (length(theta)==1) && (theta == round(theta))
            m = theta+1;
            theta = linspace(0,2*pi,m);
            theta(m) = [];
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
        hp = sctool.plotpoly(w,beta);
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
        
        % Plot circles...
        linh = gobjects(length(R),2);
        for j = 1:length(R)
            % Start with evenly spaced theta
            tp = linspace(0,2*pi,20)';
            new = true(length(tp),1);
            wp = NaN*new;
            
            if verLessThan('matlab','8.4')
                
                % The individual points will be shown as they are found
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
                neww = diskmap.dmap(zp,w,beta,z,c,qdat);
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
                [tp,wp,new] = scpadapt(tp,wp,minlen,maxlen,axlim);
                
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
            Rp = linspace(0,1,14)';
            zp = Rp*exp(1i*theta(j));
            new = true(length(zp),1);
            wp = NaN*new;
            
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
                neww = diskmap.dmap(zp(new),w,beta,z,c,qdat);
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
    
    % Used directly by hplmap/hp2disk.
    function I = dquad(z1,z2,sing1,z,beta,qdat)
        %DQUAD: Numerical quadrature for the disk map.
        %   DQUAD(ZL,ZR,S,Z,BETA,QDAT,MIDPT) performs the integration for the SC
        %   disk formula. ZL and ZR are vectors of left and right endpoints. S
        %   is a vector of integers. If ZL(k) = Z(m), then S(k) should have the
        %   value m; otherwise, S(k) should be zero.
        %
        %   Z and BETA are prevertices and turning angles for the SC map. QDAT
        %   is a matrix of quadrature data (see SCQDATA).
        %
        %   The integration is adaptive in the sense that members of Z (with
        %   nonzero BETA) that are close to the left endpoints cause
        %   subdivision. This is NOT true of singularities close to the right end.
        
        %   Copyright 1998--2001 by Toby Driscoll.
        %   $Id: dquad.m 212 2002-09-25 17:31:37Z driscoll $
        
        %   z1,z2 are vectors of left and right endpoints.  sing1 is a vector of
        %   integer indices which label the singularities in z1.  So if sing1(5)
        %   = 3, then z1(5) = z(3).  A zero means no singularity.  z is the
        %   vector of singularities; beta is the vector of associated turning
        %   angles.  qdat is quadrature data from SCQDATA.
        %
        %   Make sure that z and beta are column vectors.
        %
        %   DQUAD integrates from a possible singularity at the left end to a
        %   regular point at the right.  If both endpoints are singularities,
        %   you must break the integral into two pieces and make two calls.
        %
        %   The integral is subdivided, if necessary, so that no singularity
        %   lies closer to the left endpoint than 1/2 the length of the
        %   integration (sub)interval.
        
        nqpts = size(qdat,1);
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
            zr = za + dist*(zb-za);
            % Adjust Gauss-Jacobi nodes and weights to interval.
            ind = rem(sng+n,n+1)+1;
            nd = ((zr-za)*qdat(:,ind) + zr + za).'/2; % G-J nodes
            wt = ((zr-za)/2) * qdat(:,ind+n+1);	% G-J weights
            terms = 1 - (nd(ones(n,1),:))./bigz;
            if any(terms(:)==0)
                % Endpoints are practically coincident.
                I(k) = 0;
            else
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
                    I(k) = I(k) + exp(sum(log(1 - nd(ones(n,1),:)./bigz).*bigbeta)) * wt;
                end
            end
        end
    end
    
    % needed by CR inverse mapping
    function out = imapfun(varargin)
        out = diskmap.dimapfun(varargin{:});
    end
end

end