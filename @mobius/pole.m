function z = pole(this)

A = this.matrix;
z = double( homog( -A(2,2)/A(2,1) ) );