function bool = isinside(this,z)

bool = abs(z - this.center) < this.radius;
