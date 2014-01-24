function z = zero(this)

A = this.matrix;
z = double( homog( -A(1,2)/A(1,1) ) );