function z = point(this,t)

% Use homogeneous coordinates to define a reasonable interpolant without
% DBZ warnings
upper = 2*this.tangent*(t-1/2);
lower = 4*t.*(1-t);
z = this.base + double( homog(upper,lower) );
