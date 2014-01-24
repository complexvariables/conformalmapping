function this = zline(varargin)

this.base = [];
this.tangent = [];

switch nargin
  case 1
    if isa(varargin{1},'zline')  
      this = varargin{1}; return  % self-return
    end
    zp = varargin{1};
    % Look for two points in a vector.
    if isa(zp,'double') & length(zp)==2
      this.base = zp(1);
      this.tangent = diff(zp);
    else
      error('Expected a vector of 2 points on the line')
    end
  case 2
    % Look for a point and a tangent
    [zp,ztan] = deal(varargin{:});
    if isa(zp,'double') & isa(ztan,'double') & length(zp)==1 & length(ztan)==1
      this.base = zp;
      this.tangent = ztan;
    else
      error('Expected a point on the line and a tangent direction')
    end
  otherwise
    error('Expected a vector of 2 points or a point and a tangent')
end

this = class(this,'zline',closedcurve);

      