classdef sct
    % Class of needed and useful static routines to support S-C maps.
    
    properties
    end
    
    %% Internally necessary.
    methods (Static)
        function [wn,betan] = addVertex(w,beta,pos,window)
            %SCADDVTX Add a vertex to a polygon.
            %   [WN,BETAN] = SCADDVTX(W,BETA,POS) adds a new vertex to the polygon
            %   described by W and BETA immediately after vertex POS.  If
            %   W(POS:POS+1) are finite, the new vertex is at the midpoint of an
            %   edge; otherwise, the new vertex is a reasonable distance from its
            %   finite neighbor.
            %
            %   See also SCFIX.
            
            %   Copyright 1998 by Toby Driscoll.
            %   $Id: scaddvtx.m 298 2009-09-15 14:36:37Z driscoll $
            
            if nargin < 4, window = [-Inf Inf -Inf Inf]; end
            w = w(:);
            beta = beta(:);
            n = length(w);
            if ~pos, pos=n; end
            pos1 = rem(pos,n)+1;
            if ~any(isinf(w([pos,pos1])))	% easy case
                new = mean(w([pos,pos1]));
            else					% messy case
                % Find a pair of adjacent finite vertices as a basis.
                base = min(find(~isinf(w) & ~isinf(w([2:n,1]))));
                ang(base) = angle(w(rem(base,n)+1)-w(base));
                
                % Determine absolute angle of side pos->pos1.
                for j = [base+1:n,1:base-1]
                    ang(j) = ang(rem(j-2+n,n)+1)-pi*beta(j);
                    if j==pos, break, end
                end
                
                % Find a nice side length.
                len = abs(w([2:n,1])-w);
                avglen = mean(len(~isinf(len)));
                
                if isinf(w(pos))
                    base = w(pos1);
                    dir = exp(i*(ang(pos)+pi));
                else
                    base = w(pos);
                    dir = exp(i*(ang(pos)));
                end
                
                % Restrict point to a window (to help out graphics).
                new = base + avglen*dir;
                while real(new) < window(1) | real(new) > window(2) | ...
                        imag(new) < window(3) | imag(new) > window(4)
                    avglen = avglen / 2;
                    new = base + avglen*dir;
                end
                
            end
            
            wn = [w(1:pos);new;w(pos+1:n)];
            betan = [beta(1:pos);0;beta(pos+1:n)];
        end
        
        function beta = angles(w)
            %SCANGLE Turning angles of a polygon.
            %   SCANGLE(W) computes the turning angles of the polygon whose vertices
            %   are specified in the vector W.  The turning angle of a vertex
            %   measures how much the heading changes at that vertex from the
            %   incoming to the outgoing edge, normalized by pi.  For a finite
            %   vertex, it is equal in absolute value to (exterior angle)/pi, with a
            %   negative sign for left turns and positive for right turns.  Thus the
            %   turn at a finite vertex is in (-1,1], with 1 meaning a slit.
            %
            %   At an infinite vertex the turning angle is in the range [-3,-1] and
            %   is equal to the exterior angle of the two sides extended back from
            %   infinity, minus 2.  SCANGLE cannot determine the angle at an
            %   infinite vertex or its neighbors, and will return NaN's in those
            %   positions.
            %
            %   See also DRAWPOLY, DEMOINF.
            
            %   Copyright 1998 by Toby Driscoll.
            %   $Id: scangle.m 298 2009-09-15 14:36:37Z driscoll $
            
            w = w(:);
            n = length(w);
            if n==0
                beta = [];
                return
            end
            
            atinf = isinf(w);
            % These can't be determined.
            mask = ~(atinf | atinf([2:n,1]) | atinf([n,1:n-1]));
            
            dw = diff( w([n 1:n]) );
            dwshift = dw([2:n,1]);
            beta = NaN*ones(size(w));
            beta(mask) = angle( dw(mask).*conj(dwshift(mask)) )/pi;
            
            % It's ill-posed to tell a point (outward) from a slit (inward). Since
            % the latter is much more common and important, we'll be generous in
            % giving it the tie.
            mod = abs(beta+1) < 1e-12;
            beta(mod) = ones(size(beta(mod)));
        end
        
        function err = checkPolygon(type,w,beta,aux)
            %SCCHECK Check polygon inputs to Schwarz-Christoffel functions.
            %   SCCHECK(MAP,W,BETA) is used by the xxPARAM functions to check the
            %   validity of inputs describing the polygon to be mapped.  MAP is a
            %   string consisting of the prefix to PARAM ('d', 'dp', etc.).
            %
            %   If errors are found, execution will terminate.  Sometimes the
            %   trouble has to do with how the parameter problem is posed, which
            %   imposes a few nonobvious constraints.  The function SCFIX is
            %   provided to automatically fix such difficulties, by renumbering or
            %   perhaps adding vertices. SCCHECK output is 1 if the problem is
            %   rectifiable by SCFIX, 2 if warning status only.
            %
            %   See also SCFIX.
            
            %   Copyright 1998 by Toby Driscoll.
            %   $Id: sccheck.m 298 2009-09-15 14:36:37Z driscoll $
            
            w = w(:);
            beta = beta(:);
            n = length(w);
            atinf = logical(beta <= -1);
            renum = 1:n;
            err = 0;
            
            % Universal truths
            if length(beta)~=n
                error('Mismatched angles and vertices')
            elseif any(beta > 1) | any(beta < -3)
                error('Each angle must be in [-3,1]')
            end
            
            % Infinite vertices
            if ~strcmp(type,'de') & ~strcmp(type,'cr')
                if any(isinf(w(~atinf))) | any(~isinf(w(atinf)))
                    error('Infinite vertices must correspond to angle <= -1')
                else
                    da = diff(find(atinf));
                    if ~isempty(da) & any(da==1)
                        error('Infinite vertices must not be adjacent')
                    end
                end
            else
                if any(atinf) | any(isinf(w))
                    error('Infinite vertices not allowed')
                end
            end
            sumb = -2 * (-1)^(strcmp(type,'de'));
            
            % Orientation conventions
            if abs(sum(beta)+sumb) < 1e-9
                fprintf('\nVertices were probably given in wrong order\n')
                err = 1;
            elseif abs(sum(beta)-sumb) > 1e-9
                fprintf('\nWarning: Angles do not sum to %d\n\n',sumb)
                err = 2;
            end
            
            % Some finer points
            if strcmp(type,'hp') | strcmp(type,'d')
                if n < 3
                    error('Polygon must have at least three vertices')
                elseif any(isinf(w([1,2,n-1])))
                    fprintf('\nInfinite vertices must not be at positions 1, 2, or n-1\n')
                    err = 1;
                elseif any(abs(beta(n)-[0,1])<eps)
                    fprintf('\nSides adjacent to w(n) must not be collinear\n')
                    err = 1;
                end
            elseif strcmp(type,'cr')
                if n < 4
                    error('Polygon must have at least four vertices')
                end
            elseif strcmp(type,'de')
                if n < 2
                    error('Polygon must have at least two vertices')
                elseif (beta(n)==0 | beta(n)==1) & (n > 2)
                    fprintf('\nSides adjacent to w(n) must not be collinear\n')
                    err = 1;
                end
            elseif strcmp(type,'st')
                if n < 5
                    error('Polygon must have at least five vertices')
                end
                ends = aux;
                renum = [ends(1):n,1:ends(1)-1];
                w = w(renum);
                beta = beta(renum);
                k = find(renum==ends(2));
                if any(atinf([2,3,n]))
                    fprintf('\nVertices at (w(ends(1)) + [1,2,-1]) must be finite\n')
                    err = 1;
                elseif k-2 < 2
                    fprintf('\nThere must be at least 2 vertices between ends 1 and 2\n')
                    err = 1;
                elseif k==n
                    fprintf('\nThere must be at least one vertex between ends 2 and 1\n')
                    err = 1;
                end
            elseif strcmp(type,'r')
                corner = aux;
                renum = [corner(1):n,1:corner(1)-1];
                w = w(renum);
                beta = beta(renum);
                corner = rem(corner-corner(1)+1+n-1,n)+1;
                if n < 4
                    error('Polygon must have at least four vertices')
                elseif corner~=sort(corner)
                    error('Corners must be specified in ccw order')
                elseif isinf(w(1))
                    error('Corner(1) must be finite')
                end
                if isinf(w(2))
                    fprintf('\nVertex corner(1)+1 must be finite\n')
                    err = 1;
                end
                if any(abs(beta(n)-[0,1])<eps)
                    fprintf('\nSides adjacent to w(corner(1)-1) must not be collinear\n')
                    err = 1;
                end
                
            end
        end
        
        function [w,beta,aux] = fixPolygon(type,w,beta,aux)
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
            
            % Orientation conventions
            sumb = -2 + 4*strcmp(type,'de');
            if abs(sum(beta)+sumb) < 1e-9
                % Reverse order
                w = w([1,n:-1:2]);
                beta = scangle(w);
                renum = renum([1,n:-1:2]);
                if nargin > 3
                    aux = renum(aux);
                end
            end
            
            % Less obvious restrictions
            if strcmp(type,'hp') | strcmp(type,'d')
                shift = [2:n,1];
                % Renumber, if necessary, to meet requirements:
                %   w([1,2,n-1]) finite & sides at w(n) not collinear
                while any(isinf(w([1,2,n-1]))) | any(abs(beta(n)-[0,1])<eps)
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
                            [w,beta] = scaddvtx(w,beta,1);
                            n = n+1;
                        end
                        if isinf(w(n-1))
                            [w,beta] = scaddvtx(w,beta,n-1);
                            n = n+1;
                        end
                        renum = 1:n;
                        fprintf('\nWarning: A vertex has been added.\n\n')
                        break
                    end
                end
            elseif strcmp(type,'de')
                shift = [2:n,1];
                % Renumber, if necessary, to ensure sides at w(n) not collinear
                %   (except if n==2, which is handled explicitly anyway)
                while any(abs(beta(n)-[0,1])<eps) & (n > 2)
                    renum = renum(shift);
                    w = w(shift);
                    beta = beta(shift);
                    if renum(1)==1
                        deg = abs(beta) < eps;
                        w(deg) = [];
                        beta(deg) = [];
                        renum = 1:2;
                        n = 2;
                        fprintf('\nPolygon is a line segment; removing superfluous vertices\n\n')
                        break
                    end
                end
            elseif strcmp(type,'st')
                ends = aux;
                if isempty(ends)
                    disp('Use mouse to select images of left and right ends of the strip.')
                    figure(gcf)
                    ends = scselect(w,beta,2);
                end
                renum = [ends(1):n,1:ends(1)-1];
                w = w(renum);
                beta = beta(renum);
                k = find(renum==ends(2));
                if k < 4
                    if k < n-1
                        % Switch ends.
                        renum = [k:n,1:k-1];
                        w = w(renum);
                        beta = beta(renum);
                        k = find(renum==1);
                    else
                        % Add one or two vertices.
                        for j=1:4-k
                            [w,beta] = scaddvtx(w,beta,j);
                            n = n+1;
                            k = k+1;
                            fprintf('\nWarning: A vertex has been added.\n\n')
                        end
                    end
                end
                
                if k==n
                    % Must add a vertex in any case.
                    [w,beta] = scaddvtx(w,beta,n);
                    n = n+1;
                    fprintf('\nWarning: A vertex has been added.\n\n')
                end
                
                
                if isinf(w(2))
                    % Add two vertices.
                    for j=1:2
                        [w,beta] = scaddvtx(w,beta,j);
                        n = n+1;
                        k = k+1;
                        fprintf('\nWarning: A vertex has been added.\n\n')
                    end
                elseif isinf(w(3))
                    % Add one vertex.
                    [w,beta] = scaddvtx(w,beta,2);
                    n = n+1;
                    k = k+1;
                    fprintf('\nWarning: A vertex has been added.\n\n')
                elseif isinf(w(n))
                    [w,beta] = scaddvtx(w,beta,n);
                    n = n+1;
                    fprintf('\nWarning: A vertex has been added.\n\n')
                end
                
                aux = [1,k];
                
            elseif strcmp(type,'r')
                corner = aux;
                renum = [corner(1):n,1:corner(1)-1];
                w = w(renum);
                beta = beta(renum);
                corner = rem(corner-corner(1)+1+n-1,n)+1;
                % Note: These problems are pretty rare.
                if any(abs(beta(n)-[0,1])<eps)
                    % Try swapping sides 1-2 and 3-4.
                    if ~any(abs(beta(corner(3)-1)-[0,1])<eps) & ~isinf(w(corner(3)))
                        renum = [corner(3):n,1:corner(3)-1];
                        w = w(renum);
                        beta = beta(renum);
                        corner = sort(rem(corner-corner(3)+1+n-1,n)+1);
                    else
                        error('Collinear sides make posing problem impossible')
                    end
                end
                if isinf(w(2))
                    [w,beta] = scaddvtx(w,beta,1);
                    n = n+1;
                    corner(2:4) = corner(2:4)+1;
                end
                
                aux = corner;
            end
        end
        
        function [z0,w0] = inverseStartingPoints(prefix,wp,w,beta,z,c,qdat,aux)
            %SCIMAPZ0 (not intended for calling directly by the user)
            %   SCIMAPZ0 returns starting points for computing inverses of
            %   Schwarz-Christoffel maps.
            %
            %   Each wp(j) (in the polygon plane) requires z0(j) (in the fundamental
            %   domain) whose image w0(j) is such that the line segment from w0(j)
            %   to wp(j) lies in the target (interior or exterior) region.  The
            %   algorithm here is to choose z0(j) as a (weighted) average of
            %   successive pairs of adjacent prevertices.  The resulting w0(j) is on
            %   a polygon side.  Each choice is tested by looking for intersections
            %   of the segment with (other) sides of the polygon.
            %
            %   After randomly trying 10 weights with such prevertex pairs, the
            %   routine gives up.  Failures are pretty rare.  Slits are the most
            %   likely cause of trouble, since the intersection method doesn't know
            %   "which side" of the slit it's on.  In such a case you will have to
            %   supply starting points manually, perhaps by a continuation method.
            %
            %   See also HPINVMAP, DINVMAP, DEINVMAP, RINVMAP, STINVMAP.
            
            %   Copyright 1997 by Toby Driscoll.  Last updated 05/07/97.
            
            %   P.S. This file illustrates why the different domains in the SC
            %   Toolbox have mostly independent M-files.  The contingencies for the
            %   various geometries become rather cumbersome.
            
            n = length(w);
            shape = wp;
            wp = wp(:);
            z0 = wp;
            w0 = wp;
            from_disk = strcmp(prefix(1),'d');
            from_hp = strcmp(prefix,'hp');
            from_strip = strcmp(prefix,'st');
            from_rect = strcmp(prefix,'r');
            
            % Calculate arguments of the directed polygon edges.
            if from_strip
                % Renumber to put left end of strip first
                atinf = find(isinf(z));
                renum = [atinf(1):n 1:atinf(1)-1];
                w = w(renum);
                z = z(renum);
                beta = beta(renum);
                qdat(:,1:n) = qdat(:,renum);
                qdat(:,n+1+(1:n)) = qdat(:,n+1+renum);
                kinf = max(find(isinf(z)));
                argw = cumsum([angle(w(3)-w(2));-pi*beta([3:n,1])]);
                argw = argw([n,1:n-1]);
            else
                argw = cumsum([angle(w(2)-w(1));-pi*beta(2:n)]);
            end
            
            % Express each side in a form to allow detection of intersections.
            infty = isinf(w);
            fwd = [2:n 1];
            anchor = zeros(1,n);
            anchor(~infty) = w(~infty);
            anchor(infty) = w( fwd(infty) );        % use the finite endpoint
            direcn = exp(i*argw);
            direcn(infty) = -direcn(infty);         % reverse
            len = abs( w(fwd) - w );
            
            if from_disk
                argz = angle(z);
                argz(argz<=0) = argz(argz<=0) + 2*pi;
            end
            
            if from_rect
                % Extra argument given
                L = qdat;
                qdat = aux;
            end
            
            factor = 0.5;				% images of midpoints of preverts
            done = zeros(1,length(wp));
            m = length(wp);
            iter = 0;
            if length(qdat) > 1
                tol = 1000*10^(-size(qdat,1));
            else
                tol = qdat;
            end
            
            zbase = NaN*ones(n,1);
            wbase = NaN*ones(n,1);
            idx = [];
            while m > 0				% while some not done
                % Choose a "basis point" on each side of the polygon.
                for j = 1:n
                    if from_disk
                        if j<n
                            zbase(j) = exp(i*(factor*argz(j) + (1-factor)*argz(j+1)));
                        else
                            zbase(j) = exp(i*(factor*argz(n) + (1-factor)*(2*pi+argz(1))));
                        end
                    elseif from_hp
                        if j < n-1			% between two finite points
                            zbase(j) = z(j) + factor*(z(j+1)-z(j));
                        elseif j==n-1			% between x(n-1) & Inf
                            zbase(j) = max(10,z(n-1))/factor;
                        else				% between -Inf and x(1)
                            zbase(j) = min(-10,z(1))/factor;
                        end
                    elseif from_strip
                        if j==1
                            zbase(j) = min(-1,real(z(2)))/factor;
                        elseif j==kinf-1
                            zbase(j) = max(1,real(z(kinf-1)))/factor;
                        elseif j==kinf
                            zbase(j) = i+max(1,real(z(kinf+1)))/factor;
                        elseif j==n
                            zbase(j) = i+min(-1,real(z(n)))/factor;
                        else
                            zbase(j) = z(j) + factor*(z(j+1)-z(j));
                        end
                    elseif from_rect
                        zbase(j) = z(j) + factor*(z(rem(j,n)+1)-z(j));
                        % Can't use 0 or iK' as basis points.
                        if abs(zbase(j)) < 1e-4
                            zbase(j) = zbase(j) + .2i;
                        elseif abs(zbase(j)-i*max(imag(z))) < 1e-4
                            zbase(j) = zbase(j) - .2i;
                        end
                    end
                    
                    % Find images of all the z-plane basis points.
                    if ~from_rect
                        wbase(j) = feval([prefix,'map'],zbase(j),w,beta,z,c,qdat);
                    else
                        wbase(j) = feval([prefix,'map'],zbase(j),w,beta,z,c,L,qdat);
                    end
                    
                    % Project each base point exactly to the nearest polygon side. This
                    % is needed to make intersection detection go smoothly in borderline
                    % cases.
                    proj = real( (wbase(j)-anchor(j)) * conj(direcn(j)) );
                    wbase(j) = anchor(j) + proj*direcn(j);
                    
                end
                
                if isempty(idx)
                    % First time through, assign nearest basis point to each image point
                    [dist,idx] = min(abs( wp(~done,ones(n,1)).' - wbase(:,ones(m,1)) ));
                else
                    % Other times, just change those that failed.
                    idx(~done) = rem(idx(~done),n) + 1;
                end
                z0(~done) = zbase(idx(~done));
                w0(~done) = wbase(idx(~done));
                
                % Now, cycle thru basis points
                for j = 1:n
                    % Those points who come from basis j and need checking
                    active = (idx==j) & (~done);
                    if any(active)
                        % Assume for now that it's good
                        done(active) = ones(1,sum(active));
                        % Test line segment for intersections with other sides.
                        % We'll parameterize line segment and polygon side, compute parameters
                        % at intersection, and check parameters at intersection.
                        for k=[1:j-1,j+1:n]
                            A(:,1) = [ real(direcn(k)); imag(direcn(k)) ];
                            for p = find(active)
                                dif = (w0(p)-wp(p));
                                A(:,2) = [real(dif);imag(dif)];
                                % Get line segment and side parameters at intersection.
                                if rcond(A) < eps
                                    % All 4 involved points are collinear.
                                    wpx = real( (wp(p)-anchor(k)) / direcn(k) );
                                    w0x = real( (w0(p)-anchor(k)) / direcn(k) );
                                    if (wpx*w0x < 0) | ((wpx-len(k))*(w0x-len(k)) < 0)
                                        % Intersection interior to segment: it's no good
                                        done(p) = 0;
                                    end
                                else
                                    dif = (w0(p)-anchor(k));
                                    s = A\[real(dif);imag(dif)];
                                    % Intersection occurs interior to side? and segment?
                                    if s(1)>=0 & s(1)<=len(k)
                                        if abs(s(2)-1) < tol
                                            % Special case: wp(p) is on polygon side k
                                            z0(p) = zbase(k);
                                            w0(p) = wbase(k);
                                        elseif abs(s(2)) < tol
                                            % Happens when two sides are partially coincident (slit)
                                            % Check against normal to that side
                                            if real( conj(wp(p)-w0(p))*1i*direcn(k) ) > 0
                                                % Line segment came from "outside"
                                                done(p) = 0;
                                            end
                                        elseif s(2) > 0 & s(2) < 1
                                            % Intersection interior to segment: it's no good
                                            done(p) = 0;
                                        end
                                    end
                                end
                            end
                        end
                        
                        % Short circuit if everything is finished
                        m = sum(~done);
                        if ~m, break, end
                    end
                end
                if iter > 2*n
                    error('Can''t seem to choose starting points.  Supply them manually.')
                else
                    iter = iter + 1;
                end
                factor = rand(1);			% abandon midpoints
            end
            
            shape(:) = z0;
            z0 = shape;
            shape(:) = w0;
            w0 = shape;
        end
        
        function qdat = gaussJacobi(beta,nqpts)
            %SCQDATA Gauss-Jacobi quadrature data for SC Toolbox.
            %   SCQDATA(BETA,NQPTS) returns a matrix of quadrature data suitable for
            %   other SC routines.  BETA is a vector of turning angles corresponding
            %   to *finite* singularities (prevertices and, for exterior map, the
            %   origin).  NQPTS is the number of quadrature points per subinterval,
            %   roughly equal to -log10(error).
            %
            %   All the SC routines call this routine as needed, and the work
            %   required is small, so you probably never have to call this function
            %   directly.
            %
            %   See also GAUSSJ, HPPARAM, DPARAM, DEPARAM, STPARAM, RPARAM.
            
            %   Copyright 1998 by Toby Driscoll.
            %   $Id: scqdata.m 298 2009-09-15 14:36:37Z driscoll $
            
            n = length(beta);
            qnode = zeros(nqpts,n+1);
            qwght = zeros(nqpts,n+1);
            for j = find(beta(:)>-1)'
                [qnode(:,j),qwght(:,j)] = gaussj(nqpts,0,beta(j));
            end
            [qnode(:,n+1),qwght(:,n+1)] = gaussj(nqpts,0,0);
            qdat = [qnode,qwght];
        end
        
        function K = selectVertices(w,beta,m,titl,msg)
            %SCSELECT Select one or more vertices in a polygon.
            %   K = SCSELECT(W,BETA,M) draws the polygon given by W and BETA into
            %   the current figure window and then allows the user to select M
            %   vertices using the mouse.  If M is not given, it defaults to 1.  On
            %   exit K is a vector of indices into W.
            %
            %   SCSELECT(W,BETA,M,TITLE,MESSAGE) pops up a new figure window for the
            %   selection, with title TITLE and instructional message MESSAGE.
            %
            %   See also DRAWPOLY, PLOTPOLY, MODPOLY.
            
            %   Copyright 1998--2001 by Toby Driscoll.
            %   $Id: scselect.m 298 2009-09-15 14:36:37Z driscoll $
            
            n = length(w);
            if any(isinf(w) & isinf(w([2:n,1])))
                error('Infinite vertices must not be adjacent')
            end
            
            if nargin > 3
                fig = figure('name',titl,'numbertitle','off','integerhandle','off',...
                    'menubar','none','unit','char');
                movegui(fig,'center')
                figpos = get(fig,'pos');
                pos = [0 figpos(4)-3 figpos(3) 2.5];
                t = uicontrol('style','text','unit','char','pos',pos);
                if ~iscell(msg), msg = {msg}; end
                [str,pos2] = textwrap(t,msg);
                pos(3) = pos2(3);
                pos(1) = max(0, (figpos(3)-pos2(3))/2 );
                set(t,'pos',pos,'string',str)
                h = figpos(4)-11;
                axes('unit','char','pos',[0 4 figpos(3) h]);
            else
                fig = gcf;
            end
            
            [ehan,lhan] = plotpoly(w,beta,1);
            turn_off_hold = ~ishold;
            hold on
            
            h = lhan;
            colors = get(gca,'colororder');
            if colors(1,1) > colors(1,3)
                hilit = [0 0 1];
            else
                hilit = [1 0 0];
            end
            
            oldptr = get(gcf,'pointer');
            set(gcf,'pointer','circle');
            
            if nargin < 3
                m = 1;
            end
            
            % Begin selection(s)
            figure(fig)
            for j = 1:m
                k = [];
                while isempty(k)
                    waitforbuttonpress;
                    obj = get(gcf,'currentobj');
                    [k,tmp] = find(obj==h);
                    if isempty(k)
                        disp('Selected object not a vertex.  Try again.')
                    end
                end
                set(h(k,:),'color',hilit)
                drawnow
                K(j) = k;
            end
            set(gcf,'pointer',oldptr)
            
            % Clean up
            delete(h)
            drawnow
            if turn_off_hold
                hold off
            end
            
            if nargin > 3
                delete(fig)
            end
        end
    end  % internal methods block
    
end