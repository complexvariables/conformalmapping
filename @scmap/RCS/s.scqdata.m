h06927
s 00013/00013/00012
d D 1.2 97/04/24 09:55:08 tad 2 1
c updated documentation
c 
e
s 00025/00000/00000
d D 1.1 97/04/24 09:54:11 tad 1 0
c 
e
u
U
f b 
f e 0
t
T
I 1
function qdata = scqdata(beta,nqpts);
%SCQDATA Gauss-Jacobi quadrature data for SC Toolbox.
D 2
%       SCQDATA(BETA,NQPTS) returns a matrix of quadrature data suitable
%       for other SC routines.  BETA is a vector of turning angles
%       corresponding to *finite* singularities (prevertices and, for
%       exterior map, the origin).  NQPTS is the number of quadrature
%       points per subinterval, roughly equal to -log10(error).
%	
%	All the SC routines call this routine as needed, and the work
%	required is small, so you probably never have to call this
%	function directly.
%	
%	See also GAUSSJ, HPPARAM, DPARAM, DEPARAM, STPARAM, RPARAM.
%
%	Copyright 1996 by Toby Driscoll. Last updated 11/20/96.
E 2
I 2
%   SCQDATA(BETA,NQPTS) returns a matrix of quadrature data suitable for
%   other SC routines.  BETA is a vector of turning angles corresponding
%   to *finite* singularities (prevertices and, for exterior map, the
%   origin).  NQPTS is the number of quadrature points per subinterval,
%   roughly equal to -log10(error).
%       
%   All the SC routines call this routine as needed, and the work
%   required is small, so you probably never have to call this function
%   directly.
%       
%   See also GAUSSJ, HPPARAM, DPARAM, DEPARAM, STPARAM, RPARAM.
E 2

I 2
%   Copyright 1997 by Toby Driscoll. Last updated %G%.

E 2
n = length(beta);
qnode = zeros(nqpts,n+1);
qwght = zeros(nqpts,n+1);
for j = find(beta(:)>-1)'
  [qnode(:,j),qwght(:,j)] = gaussj(nqpts,0,beta(j));
end
[qnode(:,n+1),qwght(:,n+1)] = gaussj(nqpts,0,0);
qdata = [qnode,qwght];

E 1
