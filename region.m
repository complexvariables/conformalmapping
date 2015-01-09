classdef region < cmtobject
% REGION class.
%
% A region is defined by one or more closed curves. The interior of a
% region is defined to be to the left of the tangent vector of the closed
% curve. The exterior is the complement of the interior.
%
% A region with no boundary curves is defined to be empty in the sense of
% returning true for ISEMPTY(). This can be useful to designate a type of
% region with indeterminate geomtetry, such as a disk with unknown
% center and/or radius. It has no other consistent interpretation.
%
% r = region(p)
% r = region(p, 'interior')
%   Constructs an interior region bounded by the closed curve p or the
%   cell array of closed curves p.
% r = region(q, 'exterior')
%   Constructs an exterior region bounded by the closed curve q or the cell
%   array of closed curves q.
% r = region(p, q)
% r = region(p, 'interior', q, 'exterior')
% r = region(q, 'exterior', p, 'interior')
%   Constructs a region with p as the exterior boundary and q as the
%   interior boundary. The arguments may be cell arrays.
%
% See also closedcurve.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014,
% adapted from code by Toby Driscoll, 20??.

properties
  outerboundary
  innerboundary
end

properties(Dependent)
  m
  numinner
  numouter
end

methods
  function R = region(varargin)
    if ~nargin
      return
    end
    
    if isa(varargin{1}, 'region')
      R = varargin{1};
      return
    end
    
    switch nargin
      case 1
        p = varargin{1};
        R.outerboundary = region.checkcc(p);
        
      case 2
        [p, q] = varargin{:};
        if ischar(q)
          switch lower(q)
            case 'interior'
              R.outerboundary = region.checkcc(p);
            case 'exterior'
              R.innerboundary = region.checkcc(p);
            otherwise
              error('CMT:InvalidArgument', 'String "%s" unrecognized.', q)
          end
        else
          R.outerboundary = region.checkcc(p);
          R.innerboundary = region.checkcc(q);
        end
        
      case 4
        [p, pstr, q, qstr] = varargin{:};        
        if strcmp(pstr, 'interior') && strcmp(qstr, 'exterior')
          R.outerboundary = region.checkcc(p);
          R.innerboundary = region.checkcc(q);
        elseif strcmp(qstr, 'interior') && strcmp(pstr, 'exterior')
          R.outerboundary = region.checkcc(q);
          R.innerboundary = region.checkcc(p);
        else
          error('Invalid arguemts. Seek help.')
        end
        
      otherwise
        error('Argument configuration unrecognized.')
    end
  end % ctor
  
  function S = apply(R,f)
      ob = {};
      for i = 1:length(R.outerboundary)
          ob{i} = apply(R.outerboundary{i},f);
      end
      ib = {};
      for i = 1:length(R.innerboundary)
          ib{i} = apply(R.innerboundary{i},f);
      end
      S = region(ob,ib);
  end
  
  function b = boundary(R)
    % List region boundary as cell array of closedcurves.
    
    % Not clear the best way to do this, especially for multiply
    % connected regions. Just dump for now.
    if hasouter(R)
        b = R.outerboundary;
    else
        b = {};
    end
    if hasinner(R)
        b = [b(:); R.innerboundary(:)];
    end
    if numel(b) == 1
        b = b{1};
    end
  end
  
  function box = boundbox(R)
      zi = zeros(4*R.numinner, 1);
      for k = 1:R.numinner
          zi(4*(k - 1) + (1:4)) = cmt.bb2z(boundbox(R.innerboundary{k}));
      end
      zo = zeros(4*R.numouter, 1);
      for k = 1:R.numouter
          zo(4*(k - 1) + (1:4)) = cmt.bb2z(boundbox(R.outerboundary{k}));
      end
      box = cmt.boundbox([zi; zo]);
  end
  
  function m = connectivity(R)
      m = R.numinner + R.numouter;
  end
  
       function str = char(R)
            if isempty(R)
                str = 'empty region';
                return
            end
            
            str = 'region';
            outer = R.outerboundary;
            inner = R.innerboundary;
            
            if ~isempty(outer)
                str = [str, ' interior to:\n'];
                for k = 1:numel(outer)
                    str = [ str, '   ', char(outer{k}), '\n' ];
                end
                if ~isempty(inner)
                    str = [str, '\n and'];
                end
            end
            
            if ~isempty(inner)
                str = [str, ' exterior to:\n'];
                for k = 1:numel(inner)
                    str = [str, '   ', char(inner{k}), '\n'];
                end
            end
        end
        
        function disp(R)
            fprintf( char(R) )
            fprintf('\n')
        end

  function gd = grid(~)
    % Default: empty grid.
    warning('No grid available for this region.')
    gd = [];
  end
  
  function m = get.m(R)
      m = connectivity(R);
  end
  
  function n = get.numinner(R)
    n = numel(R.innerboundary);
  end
  
  function n = get.numouter(R)
    n = numel(R.outerboundary);
  end
  
  function fill(R, varargin)
    % Fill plot of region.
    
    newplot
    washold = ishold;
    hold on
    box on
    
    pref = cmtgetpref('graphics');
    fillargs = { pref.fillcolor,...
                'edgecolor',pref.linecolor','linewidth',pref.linewidth,...
                varargin{:} };
    
    if R.numouter > 0
        % Fill interior of the outer boundary.
        for k = 1:R.numouter
            hl = plot(R.outerboundary{k});
            fill(get(hl,'xdata'),get(hl,'ydata'),fillargs{:});
        end
        pbox = plotbox( truncate(R.outerboundary{k}) );
    else
        % Fill the plot box (as though there is no outer boundary).
        zb = zeros(4*R.numinner, 1);
        for k = 1:R.numinner
            b = truncate(R.innerboundary{k});
            zb(4*(k - 1) + (1:4)) = cmt.bb2z(plotbox(b, 1));
        end
        pbox = cmt.plotbox(zb, 2);
        fill(pbox([1 1 2 2 1]),pbox([3 4 4 3 3]),fillargs{:});
    end
    
    % Poke holes based on inner boundaries.
    if R.numinner > 0
      fillargs{1} = get(gca, 'color');
      for k = 1:R.numinner
          hl = plot(R.innerboundary{k});
          fill(get(hl,'xdata'),get(hl,'ydata'),fillargs{:});
      end
    end
    
    if ~washold
      hold off
      axis(pbox)
      aspectequal
    end
  end % fill
  
  function tf = hasgrid(~)
      % No default grid.
      tf = false;
  end
  
  function tf = hasinner(R)
    tf = ~isempty(R.innerboundary);
  end
  
  function tf = hasouter(R)
    tf = ~isempty(R.outerboundary);
  end
  
  function b = inner(R)
    if R.numinner == 1
      b = R.innerboundary{1};
    else
      b = R.innerboundary;
    end
  end
  
  function tf = isempty(R)
    % Empty region?
    tf = isempty(R.outerboundary) & isempty(R.innerboundary);
  end
  
  function tf = isexterior(R)
    % True if region has only inner boundaries. False otherwise.
    tf = hasinner(R) & ~hasouter(R);
  end
  
  function tf = isinside(R, z)
    % Is point z in region?
    if isempty(R.outerboundary)
      outin = true;
    else
      outin = isinside(R.outerboundary{1}, z);
      for k = 2:numel(R.outerboundary)
        outin = outin & isinside(R.outerboundary{k}, z);
      end
    end
    
    if isempty(R.innerboundary)
      inin = true;
    else
      inin = ~isinside(R.innerboundary{1}, z);
      for k = 2:numel(R.innerboundary)
        inin = inin & ~isinside(R.innerboundary{k}, z);
      end
    end
    
    tf = outin & inin;
  end
  
  function tf = isinterior(R)
    % True if region has only outer boundaries. False otherwise.
    tf = hasouter(R) & ~hasinner(R);
  end
  
 function tf = issimplyconnected(R)
    % True if region is simply connected (only one outer or inner boundary, but
    % not both).
    tf = (isinterior(R) & R.numouter == 1) ...
         | (isexterior(R) & R.numinner == 1);
  end
  
  function b = outer(R)
    if R.numouter == 1
      b = R.outerboundary{1};
    else
      b = R.outerboundary;
    end
  end
  
  function varargout = plot(r, varargin)
      %PLOT  Same as FILL for a region.
      [varargout{1:nargout}] = fill(r,varargin{:});

%     % Plot region without fill.
% 
%     newplot
%     washold = ishold;
%     hold on
%     
%     btag = sprintf('boundary_%s', num2hex(rand));
%     inner = r.innerboundary;
%     for k = 1:numel(inner)
%       plot(inner{k}, varargin{:}, 'tag', btag)
%     end
%     outer = r.outerboundary;
%     for k = 1:numel(outer)
%       plot(outer{k}, varargin{:}, 'tag', btag)
%     end
%     
%     if ~washold
%       hold off
%       axis(plotbox(r))
%       aspectequal
%     end
%     
%     if nargout
%       out = findobj(gca, 'tag', btag);
%     end
  end
  
  function box = plotbox(R, scale)
    if nargin < 2
      scale = [];
    end
    
    if isempty(R)
      box = [];
      return
    end
    
    box = cmt.plotbox(cmt.bb2z(boundbox(R)), scale);
  end
end

methods(Static, Hidden)
  function cc = checkcc(suitor)
    % Transform suitor to column cell of closedcurves with verification.
    if numel(suitor) == 1 && ~isa(suitor, 'cell')
      suitor = {suitor};
    end    
    for k = 1:numel(suitor)
      if ~isa(suitor{k}, 'closedcurve')
        error('Argument to region must be a closedcurve or cell array thereof.')        
      end
    end
    cc = suitor(:);
  end % checkcc
end

end
