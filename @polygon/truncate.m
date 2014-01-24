function q = truncate(p)
%TRUNCATE Truncate an unbounded polygon.
%   Q = TRUNCATE(P) returns a polygon whose finite vertices are the same
%   as those of P and whose infinite vertices have been replaced by
%   several finite ones. The new vertices are chosen by using a
%   circular "cookie cutter" on P.

%   Copyright 2002-2006 by Toby Driscoll.

w = vertex(p);
n = length(w);
tau = p.incoming;

if ~any(isinf(w)), q=p; return, end

% Find a circle containing all of the finite vertices.
wf = w(~isinf(w));
xbound = [ min(real(wf)); max(real(wf)) ];
ybound = [ min(imag(wf)); max(imag(wf)) ];
zcen = mean(xbound) + 1i*mean(ybound);
delta = norm( [ diff(xbound) diff(ybound) ]/2 );
if delta < eps, delta = 1; end
R = 2*norm(delta);

% Shift the origin to zcen.
w = w - zcen;

% Each infinite side is intersected with the circle. The infinite vertex is
% replaced by finite ones on the circle.
atinf = find(isinf(w));
v = w(1:atinf(1)-1);
for k = 1:length(atinf)
  % Indices of this, previous and next vertices. 
  m = atinf(k);
  m_prev = mod(m-2,n) + 1; 
  m_next = mod(m,n) + 1;
  % Find where the adjacent sides hit the circle.
  p1 = [ abs(tau(m))^2 2*real(tau(m)'*w(m_prev)) abs(w(m_prev))^2-R^2 ];
  t1 = roots(p1);  t1 = t1(t1>0);
  z1 = w(m_prev) + t1*tau(m);
  p2 = [ abs(tau(m_next))^2 2*real(-tau(m_next)'*w(m_next)) abs(w(m_next))^2-R^2 ];
  t2 = roots(p2);  t2 = t2(t2>0);
  z2 = w(m_next) - t2*tau(m_next);
  % Include points at intermediate angles.
  dphi = mod( angle(z2/z1), 2*pi );
  phi0 = angle(z1);
  thetanew = phi0 + unique([(0:pi/4:dphi) dphi]');
  vnew = R*exp( 1i*thetanew );
  v = [ v; vnew ];
  % Fill in finite vertices up to the next infinite one.
  if k < length(atinf)
    v = [ v; w(m+1:atinf(k+1)-1) ];
  else
    v = [ v; w(m+1:end) ];
  end
end

% Shift origin back.
v = v + zcen;
q = polygon(v);