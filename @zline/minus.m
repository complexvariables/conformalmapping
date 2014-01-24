function this = minus(this,z)

if isinf(z)
  error('Must translate by a finite number')
end

this.base = this.base - z;
