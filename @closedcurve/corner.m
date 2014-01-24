function varargout = corner(this,varargin)

v = this.corner;  % the mother struct

% Parse classes of input args.
inclass = cellfun(@class,varargin,'uniformoutput',false);

% Has an index been provided?
j = find( strcmp('double',inclass) );
if ~isempty(j)
  v = v( varargin{j} );
end

% Has a field selection been made?
j = find( strcmp('char',inclass) );
if ~isempty(j)
  varargout{1} = cat(1,v.(varargin{j}));
else
  % Output depends on nargout.
  if nargout == 3
    t = cat(1,v.param);
    z = cat(1,v.point);
    alpha = cat(1,v.alpha);
    varargout = { t,z,alpha };
  else
    varargout{1} = v;
  end
end
