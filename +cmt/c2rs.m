function [x1,x2,x3] = c2rs(z)

theta = angle(z);
absz = abs(z);
phi = atan2( absz.^2-1, 2*abs(z) );
phi(isinf(z)) = pi/2;
[x1,x2,x3] = sph2cart(theta,phi,ones(size(theta)));
