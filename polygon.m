classdef polygon < closedcurve
% POLYGON Contruct polygon object.
%   POLYGON(Z) constructs a polygon object whose vertices are the entries
%   of the complex vector Z. Each vertex should only be supplied once.
%   
%   POLYGON accepts unbounded polygons (vertices at infinity). The
%   preferred manner of doing this is by using INFVERTEX to specify the
%   directions of the edges adjacent to infinity. There must always be at
%   least one finite vertex between infinite ones.
%
%   Examples:
%
%     polygon(exp(2i*pi*(1:5)/5))  % regular pentagon
%     polygon([-1 1])  % slit (empty interior)
%     polygon( [0 infvertex(1,1i)] )  % first quadrant
%     polygon( [1i 0 infvertex(1,-1)] )  % downward step
%     polygon( [1 infvertex(1i,1i) 0 infvertex(-1i,-1i)] )  % 0 < Re z < 1
%
%   See also POLYGON/VERTEX, POLYGON/ANGLE, POLYGON/PLOT.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Adapted by Everett Kropf, 2014,
% from code by Toby Driscoll, originally 1988-20??.

properties
  vertex_
  angle_
  index_
  hvertex_
  hindex_
  incoming_
end

methods
  function P = polygon(z)
    if ~nargin
      return
    end
    
    if isa(z, 'polygon')
      P = z;
      return
    end
    
    % Find vertices and directions parallel to incoming sides.
    z = homog(z);
    m = numel(z);
    vertex = nan(size(z));
    zindex = vertex;
    incoming = vertex;
    k = 1;
    n = 0;
    lastfinite = ~isinf(z(m));
    while k <= m
      n = n + 1;
      vertex(n) = double(z(k));
      zindex(n) = k;
      if isinf(vertex(n))
        if ~lastfinite
          error('Infinite vertices cannot be adjacent.')
        end
        incoming(n) = exp(1i*angle(z(k)));
        lastfinite = false;
        k = k + 2;
      else
        k1 = mod(k - 2, m) + 1;
        if lastfinite
          incoming(n) = vertex(n) - double(z(k1));
        else
          incoming(n) = -exp(1i*angle(z(k1)));
        end
        lastfinite = true;
        k = k + 1;
      end
    end % while
    vertex = vertex(~isnan(vertex));
    zindex = zindex(~isnan(zindex));
    incoming = incoming(~isnan(incoming));
    
    % Compute the interior angles.
    atinf = isinf(vertex);
    outgoing = incoming([2:n, 1]);
    alpha = mod(angle(-incoming.*conj(outgoing))/pi, 2);
    alpha(atinf) = -mod(angle(-outgoing(atinf).*conj(incoming(atinf)))/pi, 2);
    
    % When adjacent edges are antiparallel, more testing needs ot be
    % done. In the finite case, check i fthe vertex is inside the polygon
    % defined by the others.
    isantipar = abs(alpha) < 100*eps | abs(2 - alpha) < 100*eps;
    for j = find(isantipar & ~atinf)
      if isinpoly(vertex(j), vertex([1:j-1, j+1:n]))
        alpha(j) = 2;
      else
        alpha(j) = 0;
      end
    end
    % In the inf case, truncate the infinity and check the resulting
    % triangle with its neighbors.
    for j = find(isantipar & atinf)
      jj = zindex(j);
      jp = mod(jj + 1, m) + 1;
      jm = mod(jj - 2, m) + 1;
      Z = numer(z(jj))/eps; % truncation of inf
      s = imag(-Z*conj(double(z(jp) - z(jm)) - Z));
      if s > 0
        alpha(j) = 0;
      else
        alpha(j) = -2;
      end
    end
    
    % We will always use a positive (counterclockwise) internal
    % representation.
    index = -round(sum(alpha - 1)/2);
    if index < 0
      vertex = vertex(end:-1:1);
      alpha = 2 - alpha(end:-1:1);
      index = -index;
    end
    
    n = numel(vertex);
    P.corner_ = corner({(0:n-1)'/n, vertex, alpha});
    
    P.vertex_ = vertex(:);
    P.angle_ = alpha(:);
    P.index_ = index;
    P.hvertex_ = z(:);
    P.hindex_ = zindex;
    P.incoming_ = incoming;
  end % ctor
  
  function alpha = angle(p)
    % Return polygon angles as a vector.
    
    % Original function also calculated angles, but this is done in the
    % constructor, so we're going to leave all that out. -- EK
    alpha = p.angle_;
  end
  
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
      z = p.vertex_;
      box = [min(real(z)), max(real(z)), min(imag(z)), max(imag(z))];
    else
      % We might find some finite bounds. But is there any application for this?
      box = inf*[-1 1 -1 1];
    end
  end
  
  function d = diam(p)
    %DIAM    Diameter of a polygon.
    %
    %   DIAM(P) returns max_{j,k} |P(j)-P(k)|. This may be infinite.
    
    %   Copyright 2002 by Toby Driscoll.
    %   $Id: diam.m,v 1.1 2002/09/10 19:10:41 driscoll Exp $
    
    w = p.vertex_;
    [w1, w2] = meshgrid(w);
    d = max(abs(w1(:) - w2(:)));
  end
  
  function disp(p)
    % Pretty-print a polygon.
    
    %   Copyright 1998-2003 by Toby Driscoll.
    %   $Id: display.m,v 2.4 2003/05/08 18:11:36 driscoll Exp $
    
    w = p.vertex_;
    n = numel(w);
    
    if isempty(w)
      fprintf('\n   empty polygon object\n\n')
      return
    end
    
    fprintf('\n   polygon object:\n\n')
    
    % We make disp do the heavy lifting. This way the FORMAT command works
    % here too.
    
    vstr = evalc('disp(w)');
    astr = evalc('disp(p.angle_)');
    
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
    
    fprintf(['  ' repmat('-',1,wv+4+wa) '\n']);
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
    
    if ~any(isinf(p.vertex_))
      x = p.vertex_;
      x = [real(x) imag(x)];
    else
      x = {[real(p.vertex_), imag(p.vertex_)], p.angle_};
    end
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
    
    v = p.vertex_;
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
  
  function t = isempty(p)
    %   Returns true if there are no vertices.

    %   Copyright 1998 by Toby Driscoll.
    %   $Id: isempty.m,v 2.1 1998/05/10 03:52:39 tad Exp $
    
    t = isempty(p.vertex_);
  end
  
  function tf = isinf(p)
    % Is the polygon unbounded?
    
    tf = any(isinf(p.vertex_));
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
    
    %   Copyright 1998-2003 by Toby Driscoll.
    %   $Id: isinpoly.m,v 2.3 2003/01/09 14:49:49 driscoll Exp $
    
    idx = logical(winding(p, wp, varargin{:}));
  end
  
  function in = isinside(p, z)
    % Wrapper for builtin INPOLYGON.
    
    v = p.vertex_;
    in = inpolygon(real(z), imag(z), real(v), imag(v));
  end
  
  function n = length(p)
    % Returns number of vertices, NOT the polygon boundary length.
    % FIXME: To be consistent with other boundaries, this should return
    % boundary length. Use numel(vertex(p)) instead of this function!
    
    n = numel(p.vertex_);
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
    
    %   Copyright 1998-2002 by Toby Driscoll.
    
    w = p.vertex_;
    if any(isinf(w))
      error('Invalid on unbounded polygons.')
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
  
  function r = minus(p, q)
    % Translate a polygon, or subtract the vertices of two polygons.
    r = plus(p, -q);
  end
  
  function r = mrdivide(p, q)
    % Divide a polygon by a scalar.
    if ~isa(q, 'double') || numel(q) > 1
      error('mrdivide is only defined for a scalar double.')
    end
    
    r = p;
    r.vertex_ = r.vertex_/q;
  end
  
  function r = mtimes(p, q)
    % Multiply polygon by a scalar.
    if isa(q, 'polygon')
      if isa(p, 'polygon')
        error('mtimes not defined for two polygon objects.')
      end
      [q, p] = deal(p, q);
    end
    
    r = p;
    r.vertex_ = r.vertex_*q;
  end
  
  function L = perimeter(p)
    % PERIMETER Perimeter length of a polygon.
    % (This should be the length() function. -- EK)
    
    if isinf(p)
      L = inf;
    else
      w = p.vertex_;
      L = sum(abs(diff(w([1:end, 1]))));
    end
  end
  
  function r = plus(p, q)
    % Translate a polygon, or add the vertices of two polygons.
    
    if isa(q, 'polygon')
      [q, p] = deal(p, q);
    end
    
    switch class(q)
      case 'polygon'
        if numel(q.vertex_) ~= numel(p.vertex_)
          error('Polygons mst have the same length to be added.')
        elseif isinf(p) || isinf(q)
          error('Only finite polygons may be added.')
        end
        r = polygon(p.vertex_ + q.vertex_);
        
      case 'double'
        if numel(q) > 1 && numel(q) ~= numel(p.vertex_)
          error(['Only a scalar or identical-length vector may be added ' ...
                 'to a polygon.'])
        end
        r = polygon(p.vertex_ + q(:));
    end
  end
    
  function box = plotbox(p, scale)
    % To be renamed to plotbox.
    if nargin < 2 || isempty(scale)
      scale = 1.2;
    end
    atinf = isinf(p.vertex_);
    zf = p.vertex_(~atinf);
    box = [min(real(zf)), max(real(zf)), min(imag(zf)), max(imag(zf))];
    maxdiff = max(diff(box(1:2)), diff(box(3:4)));
    if maxdiff < 100*eps
      maxdiff = 1;
    end
    fac = scale*(0.5 + 0.125*any(atinf));
    box(1:2) = mean(box(1:2)) + fac*maxdiff*[-1 1];
    box(3:4) = mean(box(3:4)) + fac*maxdiff*[-1 1];
  end

  function z = point(p, t)
    % Boundary point by parameter t in [0, 1].
    n = numel(p.vertex_);
    zc = p.vertex_;
    zh = p.hvertex_;
    zhind = p.hindex_;
    z = nan(size(t));
    t = modparam(p, t(:));
    
    % Loop by side number.
    for k = 1:n
      sidek = find(t >= (k - 1)/n & t < k/n);
      if isempty(sidek)
        continue
      end
      tk = n*t(sidek) - (k - 1);
      k1 = mod(k, n) + 1;
      if isinf(zc(k))
        tangent = numer(zh(zhind(k) + 1));
        z(sidek) = double(zc(k1) + homog((1 - tk)*tangent, tk));
      elseif isinf(zc(k1))
        tangent = numer(zh(zhind(k1)));
        z(sidek) = double(zc(k) + homog(tk*tangent, 1 - tk));
      else
        z(sidek) = interp1([0 1], zc([k k1]), tk, 'linear');
      end
    end
  end
  
  function n = size(p, m)
    % Number of vertices.
    if nargin == 1
      n = [numel(p.vertex_), 1];
    elseif m ==1
      n = numel(p.vertex_);
    else
      n = 1;
    end
  end
  
  function p = subsasgn(p, S, data)
    % Allows individual vertex assignment or property modification.
    
    if length(S) == 1 && strcmp(S.type, '()') && length(S.subs) == 1
      % Single index reference.
      p.vertex_(S.subs{1}) = data;
    else
      % Property assignment.
      if strcmp(S(1).type, '.')
        prop = S(1).subs;
        idx = strcmp([prop '_'], properties(p));
        if isempty(idx)
          error('Invalid property name %s.', prop)
        end
        if length(S) > 1 && strcmp(S(2).type, '()') && length(S(2).subs) == 1
          idx = S(2).subs{1};
        else
          idx = ':';
        end
        p.([prop '_'])(idx) = data;
      end
    end
  end
  
  function out = subsref(p, S)
    % Extract vertices by index or act as property reference.
    
    % Single index reference.
    if length(S) == 1 && strcmp(S.type, '()') && length(S.subs) == 1
      out = p.vertex_(S.subs{1});
      return
    end
    
    % Property reference.
    out = [];
    if strcmp(S(1).type, '.')
      try
        out = p.(S(1).subs);
      catch err
        if strcmp(err.identifier, 'MATLAB:noSuchMethodOrField')
          prop = S(1).subs;
          idx = strcmp([prop '_'], properties(p));
          if isempty(idx)
            error('Invalid property name %s.', prop);
          end
          out = p.([prop '_']);
        else
          rethrow(err)
        end
      end
      
      % Index on the reference.
      if length(S) > 1 && strcmp(S(2).type, '()') && length(S(2).subs) == 1
        out = out(S(2).subs{1});
      end
    end
  end
  
  function zt = tangent(p, t) %#ok<INUSD,STOUT>
    error('CMT:NotImplemented', ...
          'Placeholder function waiting on implementation.')
  end
  
  function q = truncate(p)
    % TRUNCATE Truncate an unbounded polygon.
    %   Q = TRUNCATE(P) returns a polygon whose finite vertices are the same
    %   as those of P and whose infinite vertices have been replaced by
    %   several finite ones. The new vertices are chosen by using a
    %   circular "cookie cutter" on P.
    
    %   Copyright 2002-2006 by Toby Driscoll.
    
    w = p.vertex_;
    if ~any(isinf(w))
      q = p;
      return
    end
    
    n = numel(w);
    tau = p.incoming_;
    
    % Find a circle containing all of the finite vertices.
    wf = w(~isinf(w));
    xbound = [min(real(wf)); max(real(wf))];
    ybound = [min(imag(wf)); max(imag(wf))];
    zcen = mean(xbound) + 1i*mean(ybound);
    delta = norm([diff(xbound), diff(ybound)]/2);
    if delta < eps
      delta = 1;
    end
    R = 2*norm(delta);
    
    % Shift the origin to zcen.
    w = w - zcen;
    
    % Each infinite side is intersected with the circle. The infinite vertex
    % is replaced by finite ones on the circle.
    atinf = find(isinf(w));
    v = w(1:atinf(1)-1);
    for k = 1:length(atinf)
      % Indices of this, previous, and next vertices.
      m = atinf(k);
      m_prev = mod(m - 2, n) + 1;
      m_next = mod(m, n) + 1;
      % Find where the adjacent sides hit the circle.
      p1 = [abs(tau(m))^2, 2*real(tau(m)'*w(m_prev)), abs(w(m_prev))^2 - R^2];
      t1 = roots(p1);
      t1 = t1(t1 > 0);
      z1 = w(m_prev) + t1*tau(m);
      p2 = [abs(tau(m_next))^2, 2*real(-tau(m_next)'*w(m_next)), ...
            abs(w(m_next))^2 - R^2];
      t2 = roots(p2);
      t2 = t2(t2 > 0);
      z2 = w(m_next) - t2*tau(m_next);
      % Include points at intermediate angles.
      dphi = mod(angle(z2/z1), 2*pi);
      phi0 = angle(z1);
      thetanew = phi0 + unique([(0:pi/4:dphi), dphi]');
      vnew = R*exp(1i*thetanew);
      v = [v; vnew]; %#ok<AGROW>
      % Fill in finite vertices up to the next infinite one.
      if k < length(atinf)
        v = [v; w(m + 1:atinf(k + 1) - 1)]; %#ok<AGROW>
      else
        v = [v; w(m + 1:end)]; %#ok<AGROW>
      end
    end
    
    % Shift origin back.
    v = v + zcen;
    q = polygon(v);
  end
  
  function q = uminus(p)
    %   Negate the vertices of a polygon.
    %   This may have surprising consequences if p is unbounded.
    
    %   Copyright 2003 by Toby Driscoll (driscoll@math.udel.edu).
    %   $Id: uminus.m,v 1.1 2003/03/03 16:28:04 driscoll Exp $
    
    q = polygon(-p.vertex_, p.angle_);
  end
  
  function [x, y] = vertex(p)
    %VERTEX Vertices of a polygon.
    %   VERTEX(P) returns the vertices of polygon P as a complex vector.
    %
    %   [X,Y] = VERTEX(P) returns the vertices as two real vectors.
    %
    %   See also POLYGON.
    
    %   Copyright 1998 by Toby Driscoll.
    %   $Id: vertex.m,v 2.1 1998/05/10 04:01:50 tad Exp $
    
    x = p.vertex_;
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
    
    %   Copyright 2003 by Toby Driscoll.
    %   $Id: winding.m,v 1.2 2003/01/09 14:48:25 driscoll Exp $
    
    if isinf(p)
      warning('Using a truncated version of the polygon.')
      p = truncate(p);
    end
    
    % FIXME: Which isinpoly is this supposed to be? Looks like the non-class
    % version.
    idx = double(isinpoly(wp, p.vertex_, varargin{:}));
  end
end

methods(Hidden)
  function handle = plot_(p, varargin)
    % PLOT a polygon.
    
    if isempty(p.vertex_)
      return
    end
              
    % An unbounded polygon will be truncated to make it safe for plotting.
    zplot = vertex(truncate(p));
    
    zplot = zplot([1:end 1]);
    args = plotdef.closedcurveargs;
    h = plot(real(zplot), imag(zplot), args{:});
    set(h, varargin{:});
    
    if nargout
      handle = h;
    end
  end
end

end
