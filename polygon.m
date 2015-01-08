classdef polygon < closedcurve
    %POLYGON Contruct polygon object.
    %   POLYGON(W) constructs a polygon object whose vertices are specified
    %   by the complex vector W. Cusps and cracks are allowed.
    %
    %   POLYGON(X,Y) specifies the vertices with two real vectors.
    %
    %   POLYGON(W,ALPHA) or POLYGON(X,Y,ALPHA) manually specifies the interior
    %   angles at the vertices, divided by pi.
    %
    %   POLYGON accepts unbounded polygons (vertices at infinity). However,
    %   you must supply ALPHA, and the vertices must be in counterclockwise
    %   order about the interior.
    %
    %   See also POLYGON/ANGLE, POLYGON/PLOT.
    
    % This file is a part of the CMToolkit.
    % It is licensed under the BSD 3-clause license.
    % (See LICENSE.)
    
    % Copyright Toby Driscoll, 2014.
    
    properties
        vertices = []
        angles = []
    end
    
    methods
        function P = polygon(x, y, alpha)
            if ~nargin
                return
            end
            
            if nargin < 3
                alpha = [];
            end
            
            if ~isreal(x) || nargin == 1 || (any(isinf(x)) && nargin==2)
                % Vertices passed as a complex vector
                w = x(:);
                scale = max(abs(w));
                % If first point is repeated at the end, delete the second copy
                % Thanks to Mark Embree for bug fix.
                if abs(w(end) - w(1)) < 3*eps*scale
                    w(end) = [];
                end
                if nargin > 1
                    alpha = y;
                end
            else
                % Vertices passed as two real vectors
                w = x(:) + 1i*y(:);
                scale = max(abs(w));
                % If first point is repeated at the end, delete the second copy
                if abs(w(end) - w(1)) < 3*eps*scale
                    w(end) = [];
                end
            end
            
            n = numel(w);
            if isempty(alpha)
                alpha = polygon.computeAngles(w);
            end
            index = polygon.index(alpha);
            if index < 0
                w = w([1 n:-1:2]);
                atinf = isinf(w);
                alpha = alpha([1 n:-1:2]);
                alpha(~atinf) = 2 - alpha(~atinf);
                alpha(atinf) = -2 - alpha(atinf);
            end
            
            if abs(index) > 1
                warning('CMT:BadThings', 'Polygon is multiple-sheeted.')
            end

            % Parameterize the boundary.
            function tau = tangent(t)
                thisSide = max( 1, ceil(t(:)) );
                nextSide = mod(thisSide,n) + 1;
                tau = nan(size(t));
                
                mask1 = isinf(w(nextSide));
                tau(mask1) = sign(w(nextSide(mask1)));
                mask2 = isinf(w(thisSide));
                tau(mask2) = -sign(w(thisSide(mask2)));
                mask = ~(mask1 | mask2);
                tau(mask) = w(nextSide(mask)) - w(thisSide(mask));
            end
            
            function z = position(t)
                thisSide = max( 1, ceil(t(:)) );
                nextSide = mod(thisSide,n) + 1;
                r = t(:) - thisSide + 1;   % fraction along the side
                tau = tangent(t(:));
                z = nan(size(t));
                
                % Will need the fact that infinite vertices cannot be
                % adjacent here.
                mask1 = isinf(w(nextSide));
                z(mask1) = w(thisSide(mask1)) + r(mask1)./(1-r(mask1)).*tau(mask1);
                mask2 = isinf(w(thisSide));
                z(mask2) = w(nextSide(mask2)) - (1-r(mask2))./r(mask2).*tau(mask2);
                mask = ~(mask1 | mask2);
                z(mask) = w(thisSide(mask)) + r(mask).*tau(mask);
            end
            
            P = P@closedcurve(@position,@tangent,[0 n]);
            P.vertices = w(:);
            P.angles = alpha(:);
            
        end  % constructor
        
        
        function box = boundbox(p)
            % BOUNDINGBOX Smallest box that contains the polygon.
            %   BOUNDINGBOX(P) returns the smallest box (in AXIS format) that contains
            %   the polygon P. If P is unbounded, all the entries will be infinite.
            %
            %   See also POLYGON/DIAM.
            %
            %   Copyright 2003 by Toby Driscoll.
            %   $Id: boundingbox.m,v 1.1 2003/04/25 18:46:31 driscoll Exp $
            
            if ~isinf(p)
                z = p.vertices;
                box = [min(real(z)), max(real(z)), min(imag(z)), max(imag(z))];
            else
                % We might find some finite bounds. But is there any application for this?
                box = inf*[-1 1 -1 1];
            end
        end
        
        function T = cdt(p)
            %CDT    Constrained Delaunay triangulation of polygon vertices.
            %   T = CDT(P) returns a structure representing a constrained Delaunay
            %   triangulation of the n polygon vertices. T has the fields:
            %
            %      T.edge    : 2x(2n-3) matrix of indices of edge endpoints
            %      T.triedge : 3x(n-2) matrix of triangle edge indices
            %      T.edgetri : 2x(2n-3) matrix of triangle membership indices for
            %                  the edges (boundary edges have a zero in 2nd row)
            %
            %   See also PLOTCDT.
            
            w = p.vertices;
            if any(isinf(w))
                error('CMT:NotDefined', 'CDT not possible for unbounded polygons.')
            end
            
            [e, te, et] = crtriang(w);
            [e, te, et] = crcdt(w, e, te, et);
            
            T = struct('edge', e, 'triedge', te, 'edgetri', et);
        end
        
        function d = diam(p)
            %DIAM    Diameter of a polygon.
            %
            %   DIAM(P) returns max_{j,k} |P(j)-P(k)|. This may be infinite.
            
            w = p.vertices;
            [w1, w2] = meshgrid(w);
            d = max(abs(w1(:) - w2(:)));
        end
        
        function disp(p)
            % Pretty-print a polygon.
            
            %   Copyright 1998-2003 by Toby Driscoll.
            %   $Id: display.m,v 2.4 2003/05/08 18:11:36 driscoll Exp $
            
            w = p.vertices;
            n = numel(w);
            
            if isempty(w)
                fprintf('\n   empty polygon object\n\n')
                return
            end
            
            fprintf('\n   polygon object:\n\n')
            
            % We make disp do the heavy lifting. This way the FORMAT command works
            % here too.
            
            vstr = evalc('disp(w)');
            astr = evalc('disp(p.angles)');
            
            % Parse into one cell per line.
            vc = textscan(vstr, '%s', n, 'delimiter', '\n');
            vc = vc{1};
            ac = textscan(astr, '%s', n, 'delimiter', '\n');
            ac = ac{1};
            
            % Now into matrices.
            vm = char(vc);
            am = char(ac);
            
            
            % Remove leading and trailing space blocs.
            % (Should use strtrim here? -- EK)
            idx = find(~all(vm == ' '));
            vm = vm(:,min(idx):max(idx));
            idx = find(~all(am == ' '));
            am = am(:,min(idx):max(idx));
            
            wv = max(size(vm, 2), 6);
            wa = max(size(am, 2), 8);
            b1 = blanks(2 + floor((wv - 6)/2));
            b2 = blanks(ceil((wv - 6)/2) + 4 + floor((wa - 8)/2));
            fprintf([b1 'Vertex' b2 'Angle/pi\n']);
            
            uv = min(size(vm, 2), 6);
            ua = min(size(am, 2), 8);
            b1 = blanks(2 + floor((6 - uv)/2));
            b2 = blanks(ceil((6 - uv)/2) + 4 + floor((8 - ua)/2));
            str = [repmat(b1, n, 1), vm, repmat(b2, n, 1), am];
            
            fprintf(['  ' repmat('-', 1, wv+4+wa) '\n']);
            disp(str)
            fprintf('\n\n')
        end % disp
        
        function x = double(p)
            % DOUBLE Convert polygon to double.
            %   If the polygon is bounded, DOUBLE returns the vertices in an Nx2
            %   matrix. Otherwise, it returns a cell array whose first component is
            %   the vertex matrix and whose second component is the vector of
            %   interior normalized angles.
            
            %   Copyright 1998 by Toby Driscoll.
            %   $Id: double.m,v 2.1 1998/05/10 03:51:49 tad Exp $
            
            if ~any(isinf(p.vertices))
                x = p.vertices;
                x = [real(x), imag(x)];
            else
                x = {[real(p.vertices), imag(p.vertices)], p.angles};
            end
        end
        
        function j = end(p, ~, ~)
            j = length(p);
        end
        
        function H = fill(p, varargin)
            % FILL   Plot a polygon with a filled interior.
            %   FILL(P) plots the boundary of P in blue and fills the interior of the
            %   polygon with gray. FILL(P,PROP1,VAL1,...) passes additional arguments
            %   to the built-in FILL command.
            %
            %   See also FILL.
            
            % Copyright 2003 by Toby Driscoll.
            % $Id: fill.m,v 1.3 2004/05/27 13:11:21 driscoll Exp $
            
            v = p.vertices;
            vf = v(~isinf(v));
            if any(isinf(v))
                v = vertex(truncate(p));
            end
            
            axlim = [min(real(vf)), max(real(vf)), min(imag(vf)), max(imag(vf))];
            d = max([diff(axlim(1:2)), diff(axlim(3:4))]);
            if d < eps
                d = 1;
            end
            axlim(1:2) = mean(axlim(1:2)) + 0.54*[-1 1]*d;
            axlim(3:4) = mean(axlim(3:4)) + 0.54*[-1 1]*d;
            
            % Use defaults, but allow overrides and additional settings.
            settings = {[0.75 0.75 0.85], 'edgecolor', 'b', ...
                'linewidth', 1.5, varargin{:}}; %#ok<CCAT>
            v = v([1:end, 1]);
            h = fill(real(v), imag(v), settings{:});
            
            if ~ishold
                axis equal
                axis square
                axis(axlim)
            end
            
            if nargout
                H = h;
            end
        end
        
        function [hits, loc] = intersect(p, endpt, tol)
            %INTERSECT  Find intesection of segments with polygon sides.
            %
            %   S = INTERSECT(P,ENDPT) checks for intesections between sides of the
            %   polygon P and the line segments whose endpoints are given in the
            %   complex M by 2 matrix ENDPT. If P has N sides, on return S is an
            %   M by N logical matrix with nonzeros at locations indicating
            %   intersection.
            %
            %   INTERSECT(P,ENDPT,TOL) requires that the intersection take place
            %   more than TOL away (relatively) from the segments' endpoints. By
            %   default TOL=EPS. To test truly closed segments, use
            %   INTERSECT(P,ENDPT,0); however, this is a poorly conditioned
            %   problem.
            
            n = numel(p);
            if nargin < 3
                tol = eps;
            end
            
            m = size(endpt,1);
            w = vertex(p);
            beta = angle(p)-1;
            
            % Where are the slits?
            isslit = abs(beta-1) < 2*eps;
            isslit = isslit | isslit([2:n 1]);
            
            % Find two consecutive finite vertices.
            dw = diff( w([1:n 1]) );
            K = find(~isinf(dw), 1);
            % Arguments of polygon sides.
            argw = ones(n,1);
            argw([K:n 1:K-1]) = cumsum( [angle(dw(K));-pi*beta([K+1:n 1:K-1])] );
            
            % Check each side. Solve for two parameters and check their ranges.
            hits = false(m, n);
            loc = nan(m, n);
            for k = 1:n
                tangent = exp(1i*argw(k));
                if ~isinf(w(k))
                    wk = w(k);
                    s1max = abs( w(rem(k,n)+1)-w(k) );  % parameter in [0,s1max]
                else
                    % Start from next vertex and work back.
                    wk = w(rem(k,n)+1);
                    tangent = -tangent;
                    s1max = Inf;
                end
                A(:,1) = [ real(tangent); imag(tangent) ];
                
                % Loop over the segments to be tested. The alternative is to solve a
                % block 2x2 diagonal matrix, but any collinear cases would ruin the
                % whole batch.
                for j = 1:m
                    e1e2 = endpt(j,2) - endpt(j,1);
                    A(:,2) = -[ real(e1e2); imag(e1e2) ];
                    if rcond(A) < 2*eps
                        % Segments are parallel. Check for collinearity using rotation.
                        e2 = (endpt(j,2)-wk) / tangent;
                        e1 = (endpt(j,1)-wk) / tangent;
                        if abs(imag(e1)) < 2*eps
                            % Check for overlapping.
                            x1 = min( real([e1 e2]) );
                            x2 = max( real([e1 e2]) );
                            % Do these values straddle either of the side's endpoints?
                            if (x2 >= tol) && (x1 <= s1max-tol)
                                hits(j,k) = 1;
                                loc(j,k) = wk;  % pick a place
                            end
                        end
                    else
                        % Generic case. Find intersection.
                        delta = endpt(j,1) - wk;
                        s = A \ [real(delta);imag(delta)];
                        % Check parameter ranges.
                        if s(1)>=-eps && s(1)<=s1max+eps && s(2)>=tol && s(2)<=1-tol
                            % If an end of the segment lies on a slit side, check for
                            % interior vs. exterior.
                            if isslit(k) && (abs(s(2)) < 10*eps)
                                normal = 1i*tangent;
                                if real( conj(e1e2)*normal ) < 0, break, end
                            elseif isslit(k) && (abs(s(2)-1) < 10*eps)
                                normal = 1i*tangent;
                                if real( conj(e1e2)*normal ) > 0, break, end
                            end
                            hits(j,k) = 1;
                            loc(j,k) = wk + s(1)*tangent;
                        end
                    end
                end
            end
        end
        function t = isempty(p)
            %   Returns true if there are no vertices.
            
            %   Copyright 1998 by Toby Driscoll.
            %   $Id: isempty.m,v 2.1 1998/05/10 03:52:39 tad Exp $
            
            t = isempty(p.vertices);
        end
        
        function tf = isinf(p)
            % Is the polygon unbounded?
            
            tf = any(isinf(p.vertices));
        end
        
        % FIXME: This function looks static.
        function idx = isinpoly(wp, p, varargin)
            % ISINPOLY Identify points interior/exterior to a polygon.
            %   ISINPOLY(WP,P) returns a logical vector the size of WP in which
            %   nonzero means the corresponding point is inside polygon P and zero
            %   means it is outside.
            %
            %   ISINPOLY(WP,P,TOL) considers points within TOL of the boundary to be
            %   inside P. Without this argument, points on the boundary may or may not
            %   register as inside.
            %
            %   See also POLYGON/WINDING.
            
            idx = logical(winding(p, wp, varargin{:}));
        end
        
        function in = isinside(p, z)
            % Wrapper for builtin INPOLYGON.
            
            v = p.vertices;
            in = inpolygon(real(z), imag(z), real(v), imag(v));
        end
        
        function n = length(p)
            % Returns number of vertices, NOT the polygon boundary length.
            % FIXME: To be consistent with other boundaries, this should return
            % boundary length. Use numel(vertex(p)) instead of this function!
            
            n = numel(p.vertices);
        end
        
        function [z, idx] = linspace(p, m)
            % LINSPACE Evenly spaced points around the polygon.
            %   LINSPACE(P,N) returns a vector of N points evenly spaced on the
            %   polygon P, starting with the first vertex.
            %
            %   LINSPACE(P,H) for H<1 instead uses H as an upper bound on the arc
            %   length between points.
            %
            %   [Z,IDX] = LINSPACE(...) returns the points and an identically sized
            %   vector of corresponding side indices.
            %
            %   If the polygon is unbounded, an error results.
            
            w = p.vertices;
            if any(isinf(w))
                error('CMT:NotDefined', 'Invalid on unbounded polygons.')
            end
            
            n = numel(w);
            dw = diff(w([1:n, 1]));
            
            % Arc lengths of sides.
            s = abs(dw);
            s = cumsum([0; s]);
            L = s(end);
            s = s/L; % relative arc length
            
            % Evenly spaced points in arc length.
            if m < 1
                % How many points needed?
                m = ceil(L/m) + 1;
            end
            zs = (0:m-1)'/m;
            z = zs;
            done = false(size(z));
            idx =zeros(size(z));
            
            % Translate to polygon sides.
            for j = 1:n
                mask = ~done & zs < s(j+1);
                z(mask) = w(j) + dw(j)*(zs(mask) - s(j))/(s(j+1) - s(j));
                idx(mask) = j;
                done = mask | done;
            end
        end
        
        function [p, indx] = modify(p)
            %MODIFY Modify a polygon graphically.
            %   See MODPOLY for usage instructions.
            
            [w, beta, indx] = modpoly(vertex(p), angle(p) - 1);
            p = polygon(w, beta + 1);
        end
        
        function r = minus(p, q)
            % Translate a polygon, or subtract the vertices of two polygons.
            r = plus(p, -q);
        end
        
        function r = mrdivide(p, q)
            % Divide a polygon by a scalar.
            if ~isa(q, 'double') || numel(q) > 1
                error('CMT:NotDefined', ...
                    'Right division of a polygon defined only for a scalar double.')
            end
            
            r = p;
            r.vertices = r.vertices/q;
        end
        
        function r = mtimes(p, q)
            % Multiply polygon by a scalar.
            if isa(q, 'polygon')
                if isa(p, 'polygon')
                    error('CMT:NotDefined', ...
                        'Operator "*" not defined for two polygon objects.')
                end
                [q, p] = deal(p, q);
            end
            
            r = p;
            r.vertices = r.vertices*q;
        end
        
        function L = perimeter(p)
            % PERIMETER Perimeter length of a polygon.
            
            if isinf(p)
                L = inf;
            else
                w = p.vertices;
                L = sum(abs(diff(w([1:end, 1]))));
            end
        end
        
        function h = plotcdt(p,T,varargin)
            %PLOTCDT Plot constrained Delaunay triangulation.
            %   PLOTCDT(P,T) plots the CDT of P computed by CDT. PLOTCDT(P,T,1) labels
            %   the edges and vertices.
            %
            %   H = PLOTCDT(P,T) returns a vector of handles for the edges.
            %
            %   See also CDT.
            
            han = sctool.plotptri(p.vertex, T.edge, varargin{:});
            
            if nargout > 0
                h = han;
            end
        end
        
        function box = plotbox(p, scale)
            
            if nargin < 2 || isempty(scale)
                scale = 1.2;
            end
            atinf = isinf(p.vertices);
            zf = p.vertices(~atinf);
            box = [min(real(zf)), max(real(zf)), min(imag(zf)), max(imag(zf))];
            maxdiff = max(diff(box(1:2)), diff(box(3:4)));
            if maxdiff < 100*eps
                maxdiff = 1;
            end
            fac = scale*(0.5 + 0.125*any(atinf));
            box(1:2) = mean(box(1:2)) + fac*maxdiff*[-1 1];
            box(3:4) = mean(box(3:4)) + fac*maxdiff*[-1 1];
        end
        
        function r = plus(p, q)
            % Translate a polygon, or add the vertices of two polygons.
            
            if isa(q, 'polygon')
                [q, p] = deal(p, q);
            end
            
            switch class(q)
                case 'polygon'
                    if numel(q.vertices) ~= numel(p.vertices)
                        error('Polygons mst have the same length to be added.')
                    elseif isinf(p) || isinf(q)
                        error('Only finite polygons may be added.')
                    end
                    r = polygon(p.vertices + q.vertices);
                    
                case 'double'
                    if numel(q) > 1 && numel(q) ~= numel(p.vertices)
                        error(['Only a scalar or identical-length vector may be added ' ...
                            'to a polygon.'])
                    end
                    r = polygon(p.vertices + q(:));
            end
        end
        
%         function z = point(p, t)
%             % Boundary point by parameter t in [0, 1].
%             n = numel(p.vertices);
%             zc = p.vertices;
%             zh = p.hvertex_;
%             zhind = p.hindex_;
%             z = nan(size(t));
%             t = modparam(p, t(:));
%             
%             % Loop by side number.
%             for k = 1:n
%                 sidek = find(t >= (k - 1)/n & t < k/n);
%                 if isempty(sidek)
%                     continue
%                 end
%                 tk = n*t(sidek) - (k - 1);
%                 k1 = mod(k, n) + 1;
%                 if isinf(zc(k))
%                     tangent = numer(zh(zhind(k) + 1));
%                     z(sidek) = double(zc(k1) + homog((1 - tk)*tangent, tk));
%                 elseif isinf(zc(k1))
%                     tangent = numer(zh(zhind(k1)));
%                     z(sidek) = double(zc(k) + homog(tk*tangent, 1 - tk));
%                 else
%                     z(sidek) = interp1([0 1], zc([k k1]), tk, 'linear');
%                 end
%             end
%         end
%         
        function n = size(p, m)
            % Number of vertices.
            
            if nargin == 1
                n = [numel(p.vertices), 1];
            elseif m ==1
                n = numel(p.vertices);
            else
                n = 1;
            end
        end
        
        function p = subsasgn(p, S, data)
            % Allows individual vertex assignment or property modification.
            
            if length(S) == 1 && strcmp(S.type, '()') && length(S.subs) == 1
                % Single index reference.
                p.vertices(S.subs{1}) = data;
            else
                p = builtin('subsasgn', p, S, data);
            end
        end
        
        function out = subsref(p, S)
            % Extract vertices by index or act as property reference.
            
            % Vertex reference.
            if length(S) == 1 && strcmp(S.type, '()')
                out = subsref(p.vertices, S);
                return
            end
            
            % Property reference.
            out = builtin('subsref', p, S);
        end
        
        function zt = tangent(p, t) %#ok<INUSD,STOUT>
            error('CMT:NotImplemented', ...
                'Placeholder function waiting on implementation.')
        end
        
        function [tri, x, y] = triangulate(p, h)
            %TRIANGULATE Triangulate the interior of a polygon.
            %
            %   [TRI,X,Y] = TRIANGULATE(P,H) attempts to find a reasonable
            %   triangulation of the interior of polygon P so that the typical
            %   triangle diameter is about H. If H is not specified, an automatic
            %   choice is made.
            %
            %   If P is unbounded, the polygon is first truncated to fit in a
            %   square.
            %
            %   TRIANGULATE uses DELAUNAY from Qhull, and as such does not have
            %   guaranteed success for nonconvex regions. However, things should
            %   go OK unless P has slits.
            %
            %   See also TRUNCATE, DELAUNAY.
            
            if isinf(p)
                warning('Truncating an unbounded polygon.')
                p = truncate(p);
            end
            
            w = vertex(p);
            n = length(w);
            
            if nargin < 2
                h = diam(p) / 40;
            end
            
            % Find points around boundary.
            [wb, idx] = linspace(p, h/2);   % smaller spacing to help qhull
            
            % On sides of a slit, put on extra points and perturb inward a little.
            isslit = ( abs(angle(p)-2) < 10*eps );
            slit = find( isslit | isslit([2:n 1]) );
            [wbfine, idxfine] = linspace(p, h/6);
            for k = slit(:)'
                new = (idxfine==k);
                old = (idx==k);
                wb = [ wb(~old); wbfine(new) ];  idx = [ idx(~old); idxfine(new) ];
                move = find(idx==k);
                normal = 1i*( w(rem(k,n)+1) - w(k) );
                wb(move) = wb(move) + 1e-8*normal;
            end
            
            % Take out points that are fairly close to a singularity, because it's
            % hard to find the inverse mapping there.
            for k = find(angle(p)<1)'
                close = abs(wb - w(k)) < h/3;
                wb(close) = [];  idx(close) = [];
            end
            
            % Add the polygon vertices.
            wb = [ wb; w ];
            
            % Not used? EK, 14-08-2014.
            %         idx = [ idx; (1:n)'];
            
            % Find a hex pattern that covers the interior.
            xlim = [ min(real(w)) max(real(w)) ];
            ylim = [ min(imag(w)) max(imag(w)) ];
            x = linspace(xlim(1),xlim(2),ceil(diff(xlim))/h+1);
            y = linspace(ylim(1),ylim(2),ceil(diff(ylim))/h+1);
            [X,Y] = meshgrid(x(2:end-1),y(2:end-1));
            X(2:2:end,:) = X(2:2:end,:) + (x(2)-x(1))/2;
            
            inside = isinpoly(X+1i*Y,p);
            x = [ real(wb); X(inside) ];
            y = [ imag(wb); Y(inside) ];
            
            % Triangulate using qhull.
            tri = delaunay(x,y);
            
            % Those with a boundary vertex must be examined.
            nb = length(wb);
            check = find( any( tri<=nb, 2 ) );
            
            % First, check triangle midpoints.
            idx = tri(check,:);
            z = x(idx) + 1i*y(idx);
            out = ~isinpoly( sum(z,2)/3, p );
            
            % On the rest, look for edges that cross two slit sides.
            check2 = find(~out);
            sect1 = intersect(p,z(check2,[1 2]),1e-6);
            sect2 = intersect(p,z(check2,[2 3]),1e-6);
            sect3 = intersect(p,z(check2,[3 1]),1e-6);
            out(check2( sum(sect1,2) > 1 )) = 1;
            out(check2( sum(sect2,2) > 1 )) = 1;
            out(check2( sum(sect3,2) > 1 )) = 1;
            
            tri(check(out),:) = [];
        end
        
        function q = truncate(p)
            % TRUNCATE Truncate an unbounded polygon.
            %   Q = TRUNCATE(P) returns a polygon whose finite vertices are the same
            %   as those of P and whose infinite vertices have been replaced by
            %   several finite ones. The new vertices are chosen by using a
            %   circular "cookie cutter" on P.
            
            w = p.vertices;
            n = numel(w);
            if ~any(isinf(w))
                q = p;
                return
            end
            
            % Put the finite vertices in a box.
            wf = w(~isinf(w));
            xbound = [ min(real(wf)); max(real(wf)) ];
            ybound = [ min(imag(wf)); max(imag(wf)) ];
            delta = max( diff(xbound), diff(ybound) );
            xrect = mean(xbound) + [-1;1]*delta;
            yrect = mean(ybound) + [-1;1]*delta;
            zrect = xrect([1 1 2 2]') + 1i*yrect([2 1 1 2]');
            
            % Find intersections of the box with the unbounded sides.
            [hit, loc] = intersect(p, [zrect, zrect([2:4, 1])]);
            
            % Carry over the finite vertices, inserting to substitute for the
            % infinite ones.
            atinf = find(isinf(w));
            v = w(1:atinf(1)-1);
            for k = 1:length(atinf)
                M = atinf(k);
                M1 = mod(M-2, n) + 1;
                
                % Find where the adjacent sides hit the rectangle.
                rp = loc(logical(hit(:,M1)),M1);
                %             sp = find( hit(:,M1) );
                %             rp = loc(sp,M1);
                rn = loc(logical(hit(:,M)),M);
                %             sn = find( hit(:,M) );
                %             rn = loc(sn,M);
                
                % Include the rectangle corners that are "in between".
                dt = mod( angle(rn/rp), 2*pi );
                dr = mod( angle(zrect/rp), 2*pi );
                [dr,idx] = sort(dr);
                use = dr<dt;
                v = [ v; rp; zrect(idx(use)); rn ]; %#ok<AGROW>
                if k < length(atinf)
                    v = [ v; w(M+1:atinf(k+1)-1) ]; %#ok<AGROW>
                else
                    v = [ v; w(M+1:end) ]; %#ok<AGROW>
                end
            end
            
            q = polygon(v);
        end
        
        function q = uminus(p)
            %   Negate the vertices of a polygon.
            %   This may have surprising consequences if p is unbounded.
            
            q = polygon(-p.vertices, p.angles);
        end
        
        function [x, y] = vertex(p)
            %VERTEX Vertices of a polygon.
            %   VERTEX(P) returns the vertices of polygon P as a complex vector.
            %
            %   [X,Y] = VERTEX(P) returns the vertices as two real vectors.
            %
            %   See also POLYGON.
            
            x = p.vertices;
            if nargout == 2
                y = imag(x);
                x = real(x);
            end
        end
        
        function idx = winding(p, wp, varargin)
            % WINDING Winding number of points with respect to a polygon.
            %   WINDING(P,WP) returns a vector the size of WP of winding numbers with
            %   respect to P. A zero value means the point is outside P; a value
            %   greater than 1 means it lies on multiple sheets.
            %
            %   WINDING(P,WP,TOL) makes the boundary of P "fuzzy" by a distance
            %   TOL. This may be needed to compute winding number for points on the
            %   boundary that you want to be considered "just inside."
            %
            %   See also POLYGON/ISINPOLY.
            
            if isinf(p)
                warning('CMT:BadThings', ...
                    'Using a truncated version of the polygon.')
                p = truncate(p);
            end
            
            idx = double(isinpoly(wp, p.vertices, varargin{:}));
        end
    end
    
    methods(Hidden)
        function handle = plotCurve(p, varargin)
            % PLOT a polygon.
            
            if isempty(p.vertices)
                return
            end
            
            % An unbounded polygon will be truncated to make it safe for plotting.
            zplot = vertex(truncate(p));
            
            zplot = zplot([1:end 1]);
            [cargs, pargs] = cmtplot.closedcurveArgs(varargin{:});
            h = plot(real(zplot), imag(zplot), pargs{:}, cargs{:});
            
            if nargout
                handle = h;
            end
        end
    end
    
    methods (Static)
        function alpha = computeAngles(w)
            
            % ALPHA = computeAngles(W)
            %
            % W : vector of vertices
            % ALPHA : vector of interior angles, normalized by pi
            
            n = length(w);
            
            %TODO: Fix using homog!
            if any(isinf(w))
                error('CMT:InvalidArgument', ...
                    'Cannot compute angles for unbounded polygons.')
            end
            
            % Compute angles
            incoming = w - w([n 1:n-1]);
            outgoing = incoming([2:n,1]);
            alpha = mod(angle(-incoming.*conj(outgoing))/pi ,2);
            
            % It's ill-posed to determine locally the slits (inward-pointing) from
            % points (outward-pointing). Check suspicious cases at the tips to see if
            % they are interior to the rest of the polygon.
            mask = (alpha < 100*eps) | (2-alpha < 100*eps);
            if all(mask)
                % This can happen if all vertices are collinear
                alpha(:) = 0;
                isccw = 1;				% irrelevant
                index = 1;                          % irrelevant
                return
            end
            slit = logical(isinpoly(w(mask), w(~mask)));
            fmask = find(mask);
            alpha(fmask(slit)) = 2;
            alpha(fmask(~slit)) = 0;
            
            % Now test--if incorrect, assume the orientation is clockwise
            index = sum(alpha-1)/2;
            if abs(index - round(index)) > 100*sqrt(n)*eps
                % Try reversing the interpretation of a crack
                mask = (alpha < 2*eps) | (2-alpha < 2*eps);
                alpha(~mask) = 2 - alpha(~mask);
                index = sum(alpha - 1)/2;                 
                % If still not integer, something is wrong
                if abs(index - round(index)) > 100*sqrt(n)*eps
                    error('CMT:RuntimeError', 'Invalid polygon.')
                end
            end
        end
        
        function k = index(alpha)
            % Given the angles of the polygon, determine the maximum
            % winding number about the "interior".
            k = round(sum(1-alpha)/2);  
        end
        
    end
    
end


