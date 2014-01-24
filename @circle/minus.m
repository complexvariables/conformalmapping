function this = minus(this,z)

if isinf(z)
  error('Must translate by a finite number')
end

this.center = this.center - z;
