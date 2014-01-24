function theta = angle(zeta)

% The following relies on angle(0)==0. It is standardized to [-pi,pi).
theta = mod( angle(numer(zeta)) - angle(denom(zeta)) +pi, 2*pi) - pi;