function d = dist(this,z)

% Point-circle distance

if isa(z,'circle')
  tmp = this;  this = z;  z = tmp;
end
v = z - this.center;
d = abs( abs(v) - this.radius );
