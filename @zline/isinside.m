function bool = isinside(this,z)

% Point is "inside" a line if it lies to the left as one looks along the
% tangent.
z = (z - this.base)/this.tangent;
bool = (imag(z) > 0);
