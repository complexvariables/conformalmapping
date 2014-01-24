function zr = reflect(this,z)

zr = this.center + this.radius^2/conj(z-this.center);