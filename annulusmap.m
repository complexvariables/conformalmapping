classdef annulusmap < conformalmap
    %ANNULUSMAP Conformal map to doubly connected polygonal region.
    %   ANNULUSMAP(POLYOUTER,POLYINNER) creates a Schwarz-Christoffel annulus
    %   map to the region bounded by polygons POLYOUTER and POLYINNER.
    %
    %   ANNULUSMAP(POLYOUTER,POLYINNER,'truncate') will truncate any unbounded
    %   sides. Otherwise an error is thrown for unbounded regions.
    
    %   Copyright by Toby Driscoll, 2014
    %   Written by Alfa Heryudono, 2003 and Toby Driscoll, 2014.
    
    properties
        boundary
        M, N, Z0, Z1, ALFA0, ALFA1
        isUnbounded
        qwork
        u, c, w0, w1, phi0, phi1
    end
    
    properties (Dependent)
        % backward compatibility
        ISHAPE
    end
    
    methods
        
        function map = annulusmap(target,varargin)
            
            import sctool.*
            
            % The actual domain is unknown (modulus), but we need to call
            % the parent constructor here.
            map = map@conformalmap(annulus,target);
            
            outerPolygon = outer(target);
            innerPolygon = inner(target);
            map.boundary = { outerPolygon, innerPolygon };
            
            if isinf(outerPolygon)
                if (nargin < 3) || ~isequal(varargin{1},'truncate')
                    error('Region must be bounded, or extra option given.')
                else
                    outerPolygon = truncate(outerPolygon);
                    map.isUnbounded = true;
                end
            else
                map.isUnbounded = false;
            end
            
            map.M = length(outerPolygon); map.N = length(innerPolygon);
            map.Z0 = vertex(outerPolygon).';
            map.Z1 = vertex(innerPolygon).';
            map.ALFA0 = angle(outerPolygon)';
            map.ALFA1 = 2 - angle(innerPolygon)';
            
            % Generate the Gauss-Jacobi weights & nodes
            nptq = 8;
            map.qwork = qinit(map,nptq);
            
            % Solve parameter problem.
            iguess = 0;  %(0 nonequally spaced guess or 1 equally spaced guess)
            linearc = 1; % (0 line path or 1 circular path)
            [map.u,map.c,map.w0,map.w1,map.phi0,map.phi1] = ...
                dscsolv(iguess,nptq,map.qwork,map.isUnbounded,linearc,map);
            
            % Update the map so that it has the true annulus.
            A = annulus(0,1,map.u);
            map.theDomain = A;
        end
        
        function I = get.ISHAPE(map)
            I = map.isUnbounded;
        end
        
        function disp(map)
            %DISP   Pretty-print D-SC parameters.
            %   Called automatically when a map is displayed at the command line.
            
            %   Written by Alfa Heryudono, 2003.
            
            fprintf('annulusmap:\n\n')
            if imag(map.c) < 0
                s = '-';
            else
                s = '+';
            end
            fprintf('  Conformal modulus = %5.10f        c = %.8g %c %.8gi \n\n', 1/map.u, real(map.c),s,abs(imag(map.c)));
            
            lab{1} = 'Outer';
            lab{2} = 'Inner';
            
            for k=2:-1:1;
                z1 = eval(['map.Z' int2str(k-1) '(:)']);
                alfa1 = eval(['map.ALFA' int2str(k-1) '(:)']);
                w1 = eval(['map.w' int2str(k-1) '(:)']);
                phi1 = eval(['map.phi' int2str(k-1) '(:)']);
                n = length(z1);
                
                % We make disp do the heavy lifting. This way the FORMAT command works
                % here too.
                
                z1str = evalc( 'disp(z1)' );
                alfa1str = evalc( 'disp(alfa1)' );
                w1str = evalc( 'disp(w1)' );
                phi1str = evalc( 'disp(phi1)' );
                
                % Parse into one cell per line.
                for j=1:n
                    [tmp,z1str] = strtok(z1str,sprintf('\n'));  z1c{j}=tmp;
                    [tmp,alfa1str] = strtok(alfa1str,sprintf('\n'));  alfa1c{j}=tmp;
                    [tmp,w1str] = strtok(w1str,sprintf('\n'));  w1c{j}=tmp;
                    [tmp,phi1str] = strtok(phi1str,sprintf('\n'));  phi1c{j}=tmp;
                end
                
                % Now into matrices.
                z1m = strvcat(z1c);  alfa1m = strvcat(alfa1c); w1m = strvcat(w1c);  phi1m = strvcat(phi1c);
                
                % Remove leading and trailing space blocs.
                idx = find( ~all(z1m==' ') );
                z1m = z1m(:,min(idx):max(idx));
                idx = find( ~all(alfa1m==' ') );
                alfa1m = alfa1m(:,min(idx):max(idx));
                idx = find( ~all(w1m==' ') );
                w1m = w1m(:,min(idx):max(idx));
                idx = find( ~all(phi1m==' ') );
                phi1m = phi1m(:,min(idx):max(idx));
                
                z1v = max(size(z1m,2),6);
                alfa1a = max(size(alfa1m,2),11);
                w1v = max(size(w1m,2),9);
                phi1a = max(size(phi1m,2),8);
                
                b1 = blanks(2+floor((z1v-6)/2));
                b2 = blanks(ceil((z1v-6)/2)+4+floor((alfa1a-11)/2));
                b3 = blanks(2+floor((w1v-6)/2));
                b4 = blanks(ceil((w1v-6)/2)+3+floor((phi1a-8)/2));
                b5 = blanks(floor((z1v+4+alfa1a+w1v+phi1a)/2));
                fprintf([b5 lab{k} ' Polygon' b5 '\n']);
                fprintf( [b1 'Vertex' b2 lab{k + (-1)^(k-1)} ' angle/pi' b3 'Prevertex' b4 'Angle/pi\n'] );
                %fprintf( [b1 '------' b2 '--------\n'] );
                
                uz1v = min(size(z1m,2),6);
                ualfa1a = min(size(alfa1m,2),11);
                uw1v = min(size(w1m,2),6);
                uphi1a = min(size(phi1m,2),11);
                
                b1 = blanks(2+floor((6-uz1v)/2));
                b2 = blanks(ceil((6-uz1v)/2)+4+floor((11-ualfa1a)/2));
                b3 = blanks(floor((20-uw1v)/2));
                b4 = blanks(ceil((20-uw1v)/2)+floor((5-uphi1a)/2));
                str = [ repmat(b1,n,1) z1m repmat(b2,n,1) alfa1m repmat(b3,n,1) w1m repmat(b4,n,1) phi1m];
                
                fprintf(['  ' repmat('-',1,z1v+4+alfa1a+7+w1v+4+phi1a) '\n']);
                disp(str)
                fprintf('\n\n')
                z1m = [];alfa1m = [];w1m = [];phi1m=[];
                z1c = [];alfa1c = [];w1c = [];phi1c=[];
            end
            
        end
        
        function z = eval(map,w)
            %EVAL Evaluate Schwarz-Christoffel annulus map at points.
            %   EVAL(MAP,Z) evaluates the Schwarz-Christoffel map MAP at the points Z
            %   in the canoncial annulus.
            %
            %   See also ANNULUSMAP.
            
            %   Copyright 2014 by Toby Driscoll.
            %   Written by Alfa Heryudono.
            
            % TODO: check if w is in the annulus
            
            import sctool.*

            kww = 0;
            ic = 2;
            
            % check if w is in W0.
            idx = find(map.w0 == w, 1 );
            if isempty(idx)==0
                kww = idx;
                ic = 0;
            else
                % check if w is in W1.
                idx = find(map.w1 == w, 1 );
                if isempty(idx)==0
                    kww = idx;
                    ic = 1;
                end
            end
            nptq = 8;
            
            %Making the bridge to old subroutine
            dataz = struct('M',map.M,'N',map.N,'Z0',map.Z0,'Z1',map.Z1,'ALFA0',map.ALFA0,'ALFA1',map.ALFA1,'ISHAPE',map.ISHAPE);
            z = zdsc(w,kww,ic,map.u,map.c,map.w0,map.w1,map.phi0,map.phi1,nptq,map.qwork,1,dataz);
            
        end
        function w = evalinv(map,z)
            %EVALINV   Evaluate the inverse map.
            %   EVALINV(MAP,Z) finds the inverse image of a point Z under the annulusmap
            %   MAP. That is, it maps from a doubly connected polygonal domain to the
            %   canonical annulus.
            %
            % See also ANNULUSMAP, ANNULUSMAP.EvAL.
            
            % Copyright by Toby Driscoll, 2014.
            % Written by Alfa Heryudono, 2003.
            
            % TODO: check if z is in the doubly connected region
            
            % z is not allowed to be a vertex.
            % check if z is in Z0.
            import sctool.*

            idx = find(map.Z0 == z, 1 );
            if isempty(idx)==0
                error('The point calculated is a vertex.');
            else
                % check if z is in Z1.
                idx = find(map.Z1 == z, 1 );
                if isempty(idx)==0
                    error('The point calculated is a vertex.');
                end
            end
            nptq = 8;
            eps = 1e-9;
            
            %Making the bridge to old subroutine
            dataz = struct('M',map.M,'N',map.N,'Z0',map.Z0,'Z1',map.Z1,'ALFA0',map.ALFA0,'ALFA1',map.ALFA1,'ISHAPE',map.ISHAPE);
            w = wdsc(z,map.u,map.c,map.w0,map.w1,map.phi0,map.phi1,nptq,map.qwork,eps,1,dataz);
            
            
        end
        
        function [h] = plot(f,varargin)
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
            
            if nargin==1
                type = 'curves';
            else
                type = varargin{1};
            end
            g = grid(domain(f),type);
            out = plot( apply(g,f) );
            washold = ishold;
            hold on
            
            A = range(f);
            pOuter = outer(A);
            pInner = inner(A);
            plot( pInner )
            plot( pOuter )
            axis( plotbox(pOuter,1.1))
            if ~washold
                hold off
            end
            %         end
            
            if nargout > 0
                h = out;
            end
        end
        
    end
    
    methods (Access=protected)
        % TODO: Truer vectorization?
        function z = applyMap(map,w)
            z = zeros(size(w));
            for i = 1:numel(w)
                z(i) = eval(map,w(i));
            end
        end
    end        

end
