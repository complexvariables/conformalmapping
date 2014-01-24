function zeta = cat(dim,varargin)

% Convert each argument to homog, then extract data
for n = 1:nargin-1
  h = homog(varargin{n});
  numers{n} = numer(h);  denoms{n} = denom(h);
end

% Concatenate if possible
try
  zeta = homog( cat(dim,numers{:}),cat(dim,denoms{:}) );
catch
  error('Argument dimensions are not consistent')
end