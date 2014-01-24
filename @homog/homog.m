function zeta = homog(z1,z2)

superiorto('double');

%%
% Empty matrices as default values
if nargin < 2
  z2 = [];
  if nargin < 1
    z1 = [];  
  end
end

%%
% Validate inputs
if isequal( class(z1), 'homog')
  zeta = z1;
  return
end
if ~isequal( size(z1),size(z2) )
  if isempty(z2)   % assume 1 in denominator
    z2 = ones(size(z1));
  elseif numel(z2)==1  % expand a scalar
    tmp = z1;  tmp(:) = z2;  z2 = tmp;
  else
    error('Input arguments must be scalar or of the same size')
  end
end

%%
% Transform infinities to finite representation. There is no unique choice,
% so we arbtrarily pick [\pm 1,0] based on the sign.
idx = isinf(z1);
z1(idx) = sign(z1(idx));  z2(idx) = 0;

%%
% Assign data and create object
zeta.numer = z1;
zeta.denom = z2;

zeta = class(zeta,'homog');