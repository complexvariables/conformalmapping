function s = strip(varargin)

switch nargin
  case 0
    s.tangent = [];
    s.width = [];
    s = class(s,'strip',region)
    return
  case 1
    s.tangent = varargin{1};
    s.width = 1;
  case 2
    [s.tangent,s.width] = deal(varargin{:});
  otherwise 
    error('Must give a boundary tangent, or tangent and width')
end

if ~isempty(s.tangent)
  tau = s.tangent;  z = s.width*1i*tau;
  boundary = polygon( [0 infvertex(tau,tau) z infvertex(-tau,-tau)] );
end

s = class(s,'strip',interior(boundary));
