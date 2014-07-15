classdef region
% REGION class.
%
% A region is defined by one or more closed curves. The interior of a
% region is defined to be to the left of the tangent vector of the closed
% curve. The exterior is the compliment of the interior.
%
% TBD: What does it mean for a region object to have no boundary curves?
% Currently isempty() returns true if there are no boundary curves, which
% is admittedly a bit ambiguously named. Do we mean an empty set for a
% region, or do we mean the entire plane?
%
% r = region(p)
% r = region(p, 'interiorto')
%   Constructs an interior region bounded by the closed curve p or the
%   cell array of closed curves p.
% r = region(q, 'exteriorto')
%   Constructs an exterior region bounded by the closed curve q or the cell
%   array of closed curves q.
% r = region(p, q)
% r = region(p, 'interiorto', q, 'exteriorto')
% r = region(q, 'exteriorto', p, 'interiorto')
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
  outerboundary_
  innerboundary_
end

properties(Dependent)
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
        R.outerboundary_ = region.checkcc(p);
        
      case 2
        [p, q] = varargin{:};
        if ischar(q)
          switch q
            case 'interiorto'
              R.outerboundary_ = region.checkcc(p);
            case 'exteriorto'
              R.innerboundary_ = region.checkcc(p);
            otherwise
              error('String %s unrecognized.')
          end
        else
          R.outerboundary_ = region.checkcc(p);
          R.innerboundary_ = region.checkcc(q);
        end
        
      case 4
        [p, pstr, q, qstr] = varargin{:};        
        if strcmp(pstr, 'interiorto') && strcmp(qstr, 'exteriorto')
          R.outerboundary_ = region.checkcc(p);
          R.innerboundary_ = region.checkcc(q);
        elseif strcmp(qstr, 'interiorto') && strcmp(pstr, 'exteriorto')
          R.outerboundary_ = region.checkcc(q);
          R.innerboundary_ = region.checkcc(p);
        else
          error('Invalid arguemts. Seek help.')
        end
        
      otherwise
        error('Argument configuration unrecognized.')
    end
  end % ctor
  
  function b = boundary(R)
    % List region boundary as cell array of closedcurves.
    
    % Not clear the best way to do this, especially for multiply
    % connected regions. Just dump for now.
    b = {R.outerboundary_; R.innerboundary_};
  end
  
  function disp(R)
    if isempty(R)
      fprintf('empty region\n\n')
    end
    
    fprintf('region')
    
    outer = R.outerboundary_;
    inner = R.innerboundary_;
    
    if ~isempty(outer)
      fprintf(' interior to:\n')
      for k = 1:numel(outer)
        disp(outer{k})
      end
      if ~isempty(inner)
        fprintf('\n and')
      end
    end
    
    if ~isempty(inner)
      fprintf(' exterior to:\n')
      for k = 1:numel(inner)
        disp(inner{k})
      end
    end    
  end
  
  function gd = grid(~)
    % Default empty grid.
    gd = [];
  end
  
  function n = get.numinner(R)
    n = numel(R.innerboundary_);
  end
  
  function n = get.numouter(R)
    n = numel(R.outerboundary_);
  end
  
  function fill(R, varargin)
    % Fill plot of region.
    washold = ishold;
    if ~washold
      hold on
    end
    
    fillargs = plotdef.fillargs;
    
    % Fill interiors of any outer boundaries or draw exterior region.
    if hasouter(R)
      for k = 1:R.numouter
        fill(R.outerboundary_{k}, varargin{:})
      end
    elseif isexterior(R)
      zb = zeros(4*R.numinner, 1);
      for k = 1:R.numinner
        zb(4*(k - 1) + (1:4)) = cmt.bb2z(plotbox(R.innerboundary_{k}, 1));
      end
      zb = cmt.bb2z(cmt.plotbox(zb, 2));
      fill(real(zb), imag(zb), fillargs{:}, varargin{:});
    end
    
    % Poke holes based on inner boundaries.
    if hasinner(R)
      bgcolor = get(gca, 'color');
      for k = 1:R.numinner
        fill(R.innerboundary_{k}, 'facecolor', bgcolor, ...
             'edgecolor', plotdef.filledgecolor, varargin{:})
      end
    end
            
    if ~washold
      hold off
    end
  end % fill
  
  function tf = hasinner(R)
    tf = ~isempty(R.innerboundary_);
  end
  
  function tf = hasouter(R)
    tf = ~isempty(R.outerboundary_);
  end
  
  function b = inner(R)
    if R.numinner == 1
      b = R.innerboundary_{1};
    else
      b = R.innerboundary_;
    end
  end
  
  function tf = isempty(R)
    % Empty region?
    tf = isempty(R.outerboundary_) & isempty(R.innerboundary_);
  end
  
  function tf = isexterior(R)
    % True if region has only inner boundaries. False otherwise.
    tf = hasinner(R) & ~hasouter(R);
  end
  
  function tf = isin(R, z)
    % Is point z in region?
    if isempty(R.outerboundary_)
      outin = true;
    else
      outin = isinside(R.outerboundary_{1}, z);
      for k = 2:numel(R.outerboundary_)
        outin = outin & isinside(R.outerboundary_{k}, z);
      end
    end
    
    if isempty(R.innerboundary_)
      inin = true;
    else
      inin = ~isinside(R.innerboundary_{1}, z);
      for k = 2:numel(R.innerboundary_)
        inin = inin & ~isinside(R.innerboundary_{k}, z);
      end
    end
    
    tf = outin & inin;
  end
  
  function tf = isinterior(R)
    % True if region has only outer boundaries. False otherwise.
    tf = hasouter(R) & ~hasinner(R);
  end
  
  function tf = isinside(R, z)
    % Another name for isin(R, z).
    tf = isin(R, z);
  end
  
  function tf = issimplyconnected(R)
    % True if region is simply connected (only one outer or inner boundary, but
    % not both).
    tf = (isinterior(R) & R.numouter == 1) ...
         | (isexterior(R) & R.numinner == 1);
  end
  
  function b = outer(R)
    if R.numouter == 1
      b = R.outerboundary_{1};
    else
      b = R.outerboundary_;
    end
  end
  
  function out = plot(r, varargin)
    % Plot region without fill.
    washold = ishold;
    
    newplot
    hold on
    
    btag = sprintf('boundary_%s', num2hex(rand));
    inner = r.innerboundary_;
    for k = 1:numel(inner)
      plot(inner{k}, varargin{:}, 'tag', btag)
    end
    outer = r.outerboundary_;
    for k = 1:numel(outer)
      plot(outer{k}, varargin{:}, 'tag', btag)
    end
    
    if ~washold
      hold off
    end
    
    if nargout
      out = findobj(gca, 'tag', btag);
    end
  end
  
  function box = plotbox(R, scale)
    if nargin < 2
      scale = [];
    end
    
    if isempty(R)
      box = [];
      return
    end
    
    zi = zeros(4*R.numinner, 1);
    for k = 1:R.numinner
      zi(4*(k - 1) + (1:4)) = cmt.bb2z(plotbox(R.innerboundary_{k}, 1));
    end
    zo = zeros(4*R.numouter, 1);
    for k = 1:R.numouter
      zo(4*(k - 1) + (1:4)) = cmt.bb2z(plotbox(R.outerboundary_{k}, 1));
    end
    box = cmt.plotbox([zi; zo], scale);
  end
end

methods(Static)
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
