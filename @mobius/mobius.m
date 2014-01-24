function this = mobius(varargin)
%MOBIUS Moebius transformation.
%   MOBIUS(Z,W) creates the Mobius transformation the maps the
%   3-vector Z to W. One infinity is allowed in each of Z and W.
%   
%   MOBIUS(a,b,c,d) creates the transformation
%   
%         a*z  +  b
%         ---------
%         c*z  +  d
%         
%   MOBIUS([a b; c d]) is also allowed. In either of these cases, a,b,c,d
%   should be finite complex numbers.

superiorto('double');

this.matrix = [];

switch nargin
 case 1
    A = varargin{1};
    if isa(A,'mobius')
      this = A;  return   % self-return
    end
    if isa(A,'double') & size(A)==[2 2]
      this.matrix = A;
    end
  case 4
    this.matrix = reshape( cat(1,varargin{1:4}), [2 2] ).';
  case 2
    [z,w] = deal(varargin{1:2});
    if (isa(z,'circle') || isa(z,'zline')) && ...
        (isa(w,'circle') || isa(w,'zline'))
      z = point(z,[.25 .5 .75]);
      w = point(w,[.25 .5 .75]);
    end
    if (isa(z,'double') && length(z)==3) && (isa(w,'double') && length(w)==3)
      A1 = standardmap(z);
      A2 = standardmap(w);
      this.matrix = A2\A1;
    else
      error('Expecting two vectors of 3 points or two generalized circles')
    end
end

if rcond(this.matrix) < eps
  warning('Mobius map appears to be singular')
end

this = class(this,'mobius');

end  % main function


%%
function A = standardmap(z)

  % Find the matrix of the MT from z(1:3) to [0,1,Inf]
  
  if isinf(z(1))
    A = [ 0 z(2)-z(3); 1 -z(3) ];
  elseif isinf(z(2))
    A = [ 1 -z(1); 1 -z(3) ];
  elseif isinf(z(3))
    A = [ 1 -z(1); 0 z(2)-z(1) ];
  else
    rms = z(2)-z(3);  rmq = z(2)-z(1);
    A = [ rms -z(1)*rms; rmq -z(3)*rmq ];
  end

end  % function standardmap

