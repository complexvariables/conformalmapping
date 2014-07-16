classdef (InferiorClasses = {?double}) diskmap < conformalmap
    
    properties
        polygon = []
        prevertex = []
        constant = []
    end
    
    properties (Dependent)
        conformalCenter
        accuracy
    end
    
    properties (Hidden)
        options = []
        qdata = []
    end
               
    
    %% For dependent properties
    methods
        function wc = get.conformalCenter(map)
            %CENTER Conformal center of Schwarz-Christoffel disk map.
            %   CENTER(M) returns the conformal center (image of 0) of the
            %   Schwarz-Christoffel disk map represented by M.
            %
            %   See also DISKMAP.
            wc = evaluate(map,0);
        end
            
        function  map = set.conformalCenter(map,wc)
            %   CENTER(F,WC) computes a map conformally equivalent to F but with
            %   conformal center WC (provided WC is inside the polygon of F), and
            %   returns the new map. If WC is empty, you will be asked to select it
            %   graphically.
            %
            %   See also DISKMAP.
            
            p = map.polygon;
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
            zc = dinvmap(wc,w,beta,z,map.constant,map.qdata);
                
            % Use Moebius transform to reset prevertices
            y = ((1-zc')/(1-zc))*(z-zc)./(1-zc'*z);
            y(length(y)) = 1;			% force it to be exact
            y = y./abs(y);
            
            % Recalculate constant
            mid = mean(y(1:2));
            I = diskmap.quadrature(y(1),mid,1,y,beta,map.qdata)...
                - diskmap.quadrature(y(2),mid,2,y,beta,map.qdata);
            c = diff(w(1:2))/I;
            
            % Assign new values
            map.prevertex = y;
            map.constant = c;
        end
       
        % Provided for backward compatibility. 
        function tol = get.accuracy(map)
            tol = estimateAccuracy(map);
        end
    end
    
    %% Construction.
    methods
        
        function map = diskmap(varargin)
            %DISKMAP Schwarz-Christoffel disk map object.
            %   DISKMAP(P) constructs a Schwarz-Christoffel disk map object for the
            %   polygon P. The parameter problem is solved using default options for
            %   the prevertices and the multiplicative constant.
            %
            %   DISKMAP(P,PARAM1,VALUE1,...) allows passing parameter/value
            %   pairs to affect the construction process. 
            %
            %   DISKMAP(P,Z,...) creates a diskmap object having the given prevertices Z
            %   (the mulitiplicative constant is found automatically). There is no
            %   checking to ensure that the prevertices are consistent with the
            %   given polygon. DISKMAP(P,Z,C) also uses the supplied constant. An
            %   OPTIONS argument can be added, although only the error tolerance
            %   will be used.
            %
            %   DISKMAP(Z,ALPHA,...) creates a map using the given prevertices and the
            %   interior polygon angles described by ALPHA (see POLYGON help). The
            %   image polygon is deduced by computing S-C integrals assuming a
            %   multiplicative constant of 1. DISKMAP(Z,ALPHA,C) uses the given
            %   constant instead.
            %
            %   See also DISKMAP.SET, POLYGON, CONFORMALMAP.
                                   
           z = [];
            
            % Branch based on class of first argument
            switch class(varargin{1})
                 
                case 'polygon'
                    poly = varargin{1};
                    varargin(1) = [];
                   
                case 'double'
                    % Args are the prevertex vector, then angle vector
                    z = varargin{1};
                    alpha = varargin{2};
                    varargin(1:2) = [];
                    poly = polygon(NaN*alpha*1i,alpha);
                    c = 1;
                   
                otherwise
                    msg = 'Expected "%s" to be a polygon or prevertex vector.';
                    error('CMT:diskmap:construct',...
                        msg,inputname(1))
            end 
            
            % Parse optional arguments
            for j = 1:length(varargin)
                arg = varargin{j};
                % Each arg is z or c, until options start
                if ischar(arg)
                    break
                elseif length(arg) == length(poly)
                    z = arg;
                    z = z(:);
                elseif length(arg) == 1
                    c = arg;
                else
                    msg = 'Unable to parse argument ''%s''.';
                    error('CMT:diskmap:construct',...
                        msg,inputname(j))
                end
            end
            
            % Anything left are param/value pairs. 
            % Prepare to parse them. These defaults don't matter.
            par = inputParser;
            addParamValue(par,'trace',[]);
            addParamValue(par,'tol',[]);
            addParamValue(par,'initial',[])
            
            opt = diskmap.get;   % current default options 
            parse(par,opt,varargin{j:end});           
            
            % Take actions based on what needs to be filled in
            
            if isempty(z)
                % Solve parameter problem
                % Enforce solver rules
                [w,beta] = diskmap.fixPolygon(vertex(poly),angle(poly)-1);
                poly = polygon(w);             % in case polygon was altered
                
                [z,c,qdata] = diskmap.findParameters(w,beta,...
                    opt.initial,opt.trace,opt.tol);
            end
            
            if isempty(qdata)
                % Base accuracy of quadrature on given options
                nqpts = ceil(-log10(opt.Tolerance));
                qdata = sct.gaussJacobi(angle(poly)-1,nqpts);
            end
            
            if isempty(c)
                % Find constant
                w = vertex(poly);
                beta = angle(poly)-1;
                idx = 1 + find(~isinf(z(2:end)), 1 );
                mid = mean(z([1 idx]));
                I = diskmap.quadrature(z(1),mid,1,z,beta,qdata)...
                    - diskmap.quadrature(z(idx),mid,idx,z,beta,qdata);
                c = diff(w([1 idx]))/I;
            end
            
            map = map@conformalmap(unitdisk,region(poly));
            map.polygon = poly;
            map.prevertex = z;
            map.constant = c;
            map.qdata = qdata;
            map.options = opt;
                      
            % If the polygon was not known, find it from the map.
            if any(isnan(vertex(poly)))
                map.polygon = diskmap.forwardpoly(map);
            end
                       
        end
        
        function map = continuation(map,newPolygon,varargin)
            % Continuation of given map to a new polygon
            opt = map.options;
            z0 = map.prevertex;
            if (length(z0) ~= length(newPolygon))
                error('CMT:diskmap:construct',...
                    'Polygon %s must have the same length as %s.',...
                    inputname(2),inputname(1))
            end
            if (nargin > 2)
                opt = scmapopt(opt,varargin{:});
            end
            opt = scmapopt(opt,'initial',z0);
            map = diskmap(newPolygon,z0,opt);
        end       
        

    end
    
    %% Extracting information.
    methods
        
        function acc = estimateAccuracy(map)
            %ACCURACY Apparent accuracy of a Schwarz-Christoffel disk map.
            %   ACCURACY(MAP) estimates the accuracy of the Schwarz-Christoffel disk
            %   map MAP. The technique used is to compare the differences between
            %   successive finite vertices to the integral between the corresponding
            %   prevertices, and return the maximum.
            %
            %   See also DISKMAP.
                       
            % Get data for low-level functions
            w = vertex(map.polygon);
            beta = angle(map.polygon) - 1;
            z = map.prevertex;
            c = map.constant;
            
            % Test accuracy by integrating between consecutive finite prevertices, and
            % comparing to differences of vertices.
            
            idx = find(~isinf(w));
            wf = w(idx);				% finite vertices
            
            % The two columns hold endpoint indices for integrations.
            idx = [idx(1:end) idx([2:end 1])];
            
            % Always use the center as the integration midpoint.
            mid = zeros(length(idx),1);
            
            % Do the integrations.
            I = diskmap.quadrature(z(idx(:,1)),mid,idx(:,1),z,beta,map.qdata) - ...
                diskmap.quadrature(z(idx(:,2)),mid,idx(:,2),z,beta,map.qdata);
            
            acc = max(abs( c*I - diff(wf([1:end 1])) ));
            
        end
        
        function out = char(map)
            %CHAR   Pretty-print a Schwarz-Christoffel disk map.
            
            p = map.polygon;
            w = vertex(p);
            alpha = angle(p);
            z = map.prevertex;
            c = map.constant;
            
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
            
            wc = map.conformalCenter;
            if imag(wc) < 0
                s = '-';
            else
                s = '+';
            end
            L{end+1}=sprintf('  Conformal center at %.4f %c %.4fi',real(wc),s,abs(imag(wc)));
            
            L{end+1} = sprintf('  Apparent accuracy is %.2e',map.estimateAccuracy());
            L{end+1} = ' ';
            
            out = L;
            
            
        end
        
         function wp = evaluate(map,zp,tol)
            %EVAL Evaluate Schwarz-Christoffel disk map at points.
            %   EVAL(M,ZP) evaluates the Schwarz-Christoffel map M at the points ZP
            %   in the unit disk. The default tolerance of M is used.
            %
            %   EVAL(M,ZP,TOL) attempts to give an answer accurate to TOL. If TOL is
            %   less than the accuracy of M, this is unlikely to be met.
            %
            %   See also DISKMAP, EVALINV.
            
            %   Copyright 1998 by Toby Driscoll.
            %   $Id: eval.m,v 2.1 1998/05/10 04:10:58 tad Exp $
            
            if nargin < 3
                qdata = map.qdata;
            else
                qdata = tol;
            end
            
            p = map.polygon;
            wp = NaN(size(zp));
            idx = abs(zp) <= 1+eps;  % only points in the disk
            wp(idx) = diskmap.feval(zp(idx),vertex(p),angle(p)-1,map.prevertex,map.constant,qdata);
            
        end
        
        function fp = evaluateDiff(map,zp)
            %EVALDIFF Derivative of Schwarz-Christoffel disk map at points.
            %   EVALDIFF(M,ZP) computes the derivative of the Schwarz-Christoffel
            %   disk map M at the points ZP.
            %
            %   See also DISKMAP, EVAL.
            
            %   Copyright 1998 by Toby Driscoll.
            %   $Id: evaldiff.m,v 2.1 1998/05/10 04:11:07 tad Exp $
            
            z = map.prevertex;
            c = map.constant;
            beta = angle(map.polygon) - 1;
            
            fp = diskmap.deriv(zp,z,beta,c);
        end
        
        function [zp,flag] = evaluateInv(map,wp,tol,z0)
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
            
            %   Copyright 1998 by Toby Driscoll.
            %   $Id: evalinv.m,v 2.2 2002/09/05 20:06:41 driscoll Exp $
            
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
                qdata = map.qdata;
                tol = estimateAccuracy(map);
            end
            
            if ~isempty(z0)
                if length(z0) == 1
                    z0 = repmat(z0,size(wp));
                elseif any(size(z0) ~= size(wp))
                    msg = 'Argument %s must be a complex scalar or the same size as %s.';
                    error(sprintf(msg,inputname(z0),inputname(1)));
                end
            end
            
            p = map.polygon;
            [zp,flag] = dinvmap(wp,vertex(p),angle(p)-1,map.prevertex,map.constant,...
                qdata,z0,[0 tol]);
        end
        
    end       
      
    %% Static.
    methods (Static)
        
        function p = forwardpoly(map)
            %   Given a diskmap M, FORWARDPOLY(M) returns the polygon that is
            %   formed using the prevertices, angles, and quadrature data of that
            %   map. If the prevertices were found from the solution of a
            %   parameter problem, then the result should agree closely with the
            %   original polygon that was supplied.
            
            %   Copyright (c) 1998 by Toby Driscoll.
            %   $Id: forwardpoly.m,v 1.1 1998/06/22 22:32:30 tad Exp $
            
            z = map.prevertex;
            alpha = angle(polygon(map));
            c = map.constant;
            
            n = length(z);
            
            % Since there is no parameter problem, use high accuracy in quadrature.
            qdata = sct.gaussJacobi(alpha-1,16);
            
            w = zeros(n,1);
            atinf = (alpha < eps);
            w(atinf) = Inf;
            
            % Endpoints of integrations
            idx = find(~atinf);
            idx = [idx(1:end-1) idx(2:end)];
            
            % Origin is midpoint of every integration
            mid = zeros(length(idx),1);
            
            % Integrations
            I = diskmap.quadrature(z(idx(:,1)),mid,idx(:,1),z,alpha-1,qdata) - ...
                diskmap.quadrature(z(idx(:,2)),mid,idx(:,2),z,alpha-1,qdata);
            
            % Deduce vertices
            w(~atinf) = c*cumsum([0;I]);
            
            p = polygon(w,alpha);
            
        end
        
        function prefs = get(varargin)
            prefs = cmt.get('diskmap',varargin{:});
            if isempty(prefs)
                % Set up the factory defaults.
                diskmap.set;
                prefs = cmt.get('diskmap',varargin{:});
            end
        end
        
        function set(varargin)
            %
            if (nargin == 0) || isempty(cmt.get('diskmap'))
                % Factory defaults.
                varargin = { ...
                    'trace','partial', ...
                    'tol',1e-12, ...
                    'initial',[],...
                    varargin{:} };
            end
            
            cmt.set('diskmap',varargin{:})  
        end
 
        
    end
    
    %% Private.
    methods (Static,Access=private)
        
        function fprime = deriv(zp,z,beta,c)
            %DDERIV Derivative of the disk map.
            %   DDERIV(ZP,Z,BETA,C) returns the derivative at the points of ZP of
            %   the Schwarz-Christoffel disk map defined by Z, BETA, and C.
            %
            %   See also DPARAM, DMAP.
            
            %   Copyright 1998 by Toby Driscoll.
            %   $Id: dderiv.m 7 1998-05-10 04:37:19Z tad $
            
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
        function wp = feval(zp,w,beta,z,c,qdat)
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
            
            %   Copyright 1998 by Toby Driscoll.
            %   $Id: dmap.m 7 1998-05-10 04:37:19Z tad $
            
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
                qdat = cmt.gaussJacobi(beta,8);
            elseif length(qdat)==1
                tol = qdat;
                qdat = cmt.gaussJacobi(beta,max(ceil(-log10(tol)),8));
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
                    wc = w(n-1) + c*diskmap.quadrature(z(n-1),0,n-1,z,beta,qdat);
                else
                    wc = w(n) + c*diskmap.quadrature(z(n),0,n,z,beta,qdat);
                end
            end
            
            % zs = the starting singularities
            zs = z(sing);
            % ws = map(zs)
            ws = w(sing);
            
            % Compute the map directly at "normal" points.
            normal = ~bad & ~vertex;
            if any(normal)
                I = diskmap.quadrature(zs(normal),zp(normal),sing(normal),z,beta,qdat);
                wp(normal) = ws(normal) + c*I;
            end
            
            % Compute map at "bad" points, using conformal center as basis, to avoid
            % integration where right endpoint is too close to a singularity.
            if any(bad)
                I = diskmap.quadrature(zp(bad),zeros(sum(bad),1),zeros(sum(bad),1),z,beta,qdat);
                wp(bad) = wc - c*I;
            end
            
            wp = reshape(wp,shape);
        end
        function [z,c,qdat] = findParameters(w,beta,z0,trace,tol)
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
            
            import nesolver.nesolve
            n = length(w);				% no. of vertices
            w = w(:);
            beta = beta(:);
            
            method = 2;  % trust region only
            nqpts = max(ceil(-log10(tol)),4);
            qdat = sct.gaussJacobi(beta,nqpts); 		% quadrature data
            
            atinf = (beta <= -1);
            
            if n==3
                % Trivial solution
                z = [-1i; (1-1i)/sqrt(2); 1];
                
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
                    Am = moebius(z0(n-2:n),[-1;-1i;1]);
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
                trace = find(strcmp(trace,{'off','partial','full'})) - 1;
                opt = zeros(16,1);
                opt(1) = trace;
                opt(2) = method;
                opt(6) = 100*(n-3);
                opt(8) = tol;
                opt(9) = min(eps^(2/3),tol/10);
                opt(12) = nqpts;
                try
                    %opt = [100*(n-3),-1,0.5,1];
                    %[y,itHist] = nsold(y0,@(y) dpfun(y,fdat),[tol,tol],opt);
                    [y,termcode] = nesolve('dpfun',y0,opt,fdat);
                catch
                    % Have to delete the "waitbar" figure if interrupted
                    close(findobj(allchild(0),'flat','Tag','TMWWaitbar'));
                    error(lasterr)
                end
                %if termcode~=1
                %    warning('Nonlinear equations solver did not terminate normally.')
                %end
                
                % Convert y values to z
                cs = cumsum(cumprod([1;exp(-y)]));
                theta = pi*cs(1:n-3)/cs(n-2);
                z = ones(n,1);
                z(1:n-3) = exp(1i*theta);
                z(n-2:n-1) = [-1;-1i];
            end
            
            % Determine scaling constant
            mid = (z(1)+z(2))/2;
            c = (w(1) - w(2))/...
                (diskmap.quadrature(z(2),mid,2,z,beta,qdat) - diskmap.quadrature(z(1),mid,1,z,beta,qdat));       
        end
        
        function [w,beta,aux] = fixPolygon(w,beta,aux)
            %SCFIX  Fix polygon to meet Schwarz-Christoffel toolbox constraints.
            %   [W,BETA] = SCFIX(MAP,W,BETA) attempts to fix a problem in the given
            %   polygon that arises from the posing of the parameter problem. SCFIX
            %   is used when a call to xxPARAM results in an error and so
            %   advises. In this case the polygon as given violates some fairly
            %   arbitrary constraint. SCFIX remedies the situation by renumbering
            %   the vertices, or, if necessary, adding a trivial (zero-turn)
            %   vertex. MAP is one of {'hp','d','de','st','r'}. If one additional
            %   input and output argument is given, it represents the indices of the
            %   strip ends or the rectangle corners.
            %
            %   See also SCCHECK, SCADDVTX.
            
            %   Copyright 1998 by Toby Driscoll.
            %   $Id: scfix.m 298 2009-09-15 14:36:37Z driscoll $
            
            %   You may wonder, why not let the xxPARAM functions call SCFIX
            %   automatically?  The trouble with that approach is that since a
            %   function can't modify its inputs in the calling workspace,
            %   either the xxPARAM functions would have to return more
            %   arguments, or the mapping and plotting functions also would have
            %   to detect and correct the problem every time they're called.
            %   The problem is rare enough that this method seems adequate.
            
            w = w(:);
            beta = beta(:);
            n = length(w);
            renum = 1:n;
            
            % Less obvious restrictions
            shift = [2:n,1];
            % Renumber, if necessary, to meet requirements:
            %   w([1,2,n-1]) finite & sides at w(n) not collinear
            while any(isinf(w([1,2,n-1]))) || any(abs(beta(n)-[0,1])<eps)
                renum = renum(shift);
                w = w(shift);
                beta = beta(shift);
                if renum(1)==1 			% tried all orderings
                    % First, be sure beta(n) is no longer a problem.
                    if all((abs(beta-1)<eps)|(abs(beta)<eps))
                        error('Polygon has empty interior!')
                    end
                    while any(abs(beta(n)-[0,1])<eps)
                        w = w(shift);
                        beta = beta(shift);
                    end
                    % Next, add one or two vertices as needed.
                    if any(isinf(w(1:2)))
                        [w,beta] = sct.addVertex(w,beta,1);
                        n = n+1;
                    end
                    if isinf(w(n-1))
                        [w,beta] = sct.addVertex(w,beta,n-1);
                        n = n+1;
                    end
                    renum = 1:n;
                    fprintf('\nWarning: A vertex has been added.\n\n')
                    break
                end
            end
            
        end
                        
        function I = quadrature(z1,z2,sing1,z,beta,qdat)
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
        
   end
    
    
end