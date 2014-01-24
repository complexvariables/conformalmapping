function this = mtimes(this,z)

% Swap for correct ordering
if isa(z,'zline')
  tmp = this;  this = z;  z = tmp;
end

if isinf(z)
  error('Must scale by a finite number')
end

this.base = this.base*z;
this.tangent = this.tangent*z;
