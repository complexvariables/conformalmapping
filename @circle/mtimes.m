function this = mtimes(this,z)

% Swap for correct ordering
if isa(z,'circle')
  tmp = this;  this = z;  z = tmp;
end

if isinf(z)
  error('Must scale by a finite number')
end

this.center = this.center*z;
this.radius = this.radius*abs(z);
