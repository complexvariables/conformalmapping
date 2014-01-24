function f = conformalmap(varargin)

f.function = [];
f.domain = [];
f.image = [];

% Empty case first
if nargin==0
  f = class(f,'conformalmap');
  return
end

% Self-return case
if nargin==1 & isa(varargin{1},'conformalmap')
  f = varargin{1};
  return
end

% The following checks input classes. 
[s{1:nargin}] = deal('conformalmap');
isconformal = cellfun(@isa,varargin,s);
[s{1:nargin}] = deal('function_handle');
isfunh = cellfun(@isa,varargin,s);
[s{1:nargin}] = deal('mobius');
ismob = cellfun(@isa,varargin,s);
[s{1:nargin}] = deal('region');
isregion = cellfun(@isa,varargin,s);

% If all inputs are themselves maps, compose them. 
if all( isconformal )
  f.function = varargin;
  f.domain = domain(varargin{1});
  f.image = image(varargin{end});
% Look for an explicitly given function in the first argument, or exactly
% three arguments.
elseif (isfunh(1) | ismob(1)) | nargin==3
  [f.function,f.domain,f.image] = deal(varargin{:});
% A domain and image only mean to deduce the map.  
elseif nargin==2 & all( isregion )
  f = findthemap(varargin{:});
  return
else
  error('Syntax does not follow one of the allowed formats')
end

f = class(f,'conformalmap');

end  % conformalmap()


function f = findthemap(domain,image)

if isa(domain,'disk')
  if isa( boundary(image), 'polygon' )
    f = diskmap(boundary(image));
  end
elseif isa(domain,'annulus')
  f = annulusmap(image);
end

end  % findthemap()