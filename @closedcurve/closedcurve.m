function cc = closedcurve(varargin)

cc.point = [];
cc.corner.param = [];
cc.corner.point = [];
cc.corner.alpha = [];

% Self-return case, or find the parameterization function.
if nargin > 0
  if isa(varargin{1},'closedcurve')
    cc = varargin{1};  return
  else
    cc.point = varargin{1};
  end
end

% Corner information
if nargin > 1
  v = varargin{2};
  % Cell and struct forms are allowed
  if isa(v,'struct')
    if ~all( isfield(v,{'param','point','alpha'}) )
      error('Corner information was supplied incorrectly')
    end
    cc.corner = orderfields(v,cc.corner);
  elseif isa(v,'cell')
    cc.corner = struct( 'param',num2cell(v{1}), ...
                        'point',num2cell(v{2}), ...
                        'alpha',num2cell(v{3}) );
  else
    error('Corner information was supplied incorrectly')
  end
end
  
cc = class(cc,'closedcurve');
