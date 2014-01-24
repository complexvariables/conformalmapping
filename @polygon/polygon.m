function p = polygon(z)
%POLYGON Contruct polygon object.
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

%   Copyright (c) 1998 by Toby Driscoll.


superiorto('double');

switch nargin
  case 0
    vertex = [];  alpha = [];  index = [];  z = homog([]);  zindex = [];
    incoming = [];
  case 1
    if isa(z,'polygon'), p = z;  return,  end
    % Find vertices and directions parallel to incoming sides
    z = homog(z);
    m = length(z(:));
    k = 1;  n = 0;  lastfinite = ~isinf(z(m));
    while k <= m
      n = n+1;
      vertex(n) = double(z(k));
      zindex(n) = k;
      if isinf(vertex(n))
        if ~lastfinite
          error('Infinite vertices cannot be adjacent')
        end
        incoming(n) = exp(1i*angle(z(k)));
        lastfinite = false;
        k = k+2;
      else
        k1 = mod(k-2,m) + 1;
        if lastfinite
          incoming(n) = vertex(n) - double(z(k1));
        else
          incoming(n) = -exp(1i*angle(z(k1)));
        end
        lastfinite = true;
        k = k+1;
      end
    end
    
    % Compute the interior angles
    atinf = isinf(vertex);
    outgoing = incoming([2:n,1]);
    alpha = mod( angle(-incoming.*conj(outgoing))/pi, 2);
    alpha(atinf) = -mod( angle(-outgoing(atinf).*conj(incoming(atinf)))/pi, 2);
    
    % When adjacent edges are antiparallel, more testing needs to be done
    % In finite case, check if the vertex is inside the polygon defined by
    % the others
    for j = find( ((abs(alpha) < 100*eps) | (abs(2-alpha) < 100*eps)) & ~atinf )
      if isinpoly( vertex(j), vertex([1:j-1 j+1:n]) );
        alpha(j) = 2;
      else
        alpha(j) = 0;
      end
    end
    % In the Inf case, truncate the infinity and check the resulting
    % triangle with its neighbors
    for j = find( ((abs(alpha) < 100*eps) | (abs(-2-alpha) < 100*eps)) & atinf )
      jj = zindex(j);
      jp = mod(jj+1,m) + 1;
      jm = mod(jj-2,m) + 1;
      Z = numer(z(jj))/eps;  % truncation of Inf
      s = imag( -Z*conj(double(z(jp)-z(jm))-Z) );
      if s>0
        alpha(j) = 0;
      else
        alpha(j) = -2;
      end
    end
    
    % We will always use a positive (counterclockwise) internal
    % representation
    index = -round(sum(alpha-1)/2);
    if index < 0
      vertex = vertex(end:-1:1);
      alpha = 2-alpha(end:-1:1);
      index = -index;
    end
    
end

p.vertex = vertex(:);
p.angle = alpha(:);
p.index = index;
p.hvertex = z(:);
p.hindex = zindex;
p.incoming = incoming;

n = length(p.vertex);
cc = closedcurve(@point,{(0:n-1)'/n, p.vertex, p.angle});
p = class(p,'polygon',cc);

end

