function z = double(zeta)

status = warning('query','MATLAB:divideByZero');
warning('off','MATLAB:divideByZero');

z = numer(zeta)./denom(zeta);
% Result of complex/0 is not reliably just Inf, hence:
z(isinf(z)) = Inf;

warning(status.state,'MATLAB:divideByZero');