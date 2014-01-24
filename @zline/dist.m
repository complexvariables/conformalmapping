function d = dist(this,z)

% Point-line distance

if isa(z,'zline')
  tmp = this;  this = z;  z = tmp;
end
v = z - this.base;
s = sign(1i*this.tangent);
d = abs( real(v)*real(s) + imag(v)*imag(s) );
