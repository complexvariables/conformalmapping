classdef conformalmap
% CONFORMALMAP base class.
%
% This class has two purposes. As a base class, which would normally be
% abstract, but we need its constructor to work for the second purpose,
% which is to act as an automated map creator.
%
% For example, with certain combinations of domains and ranges,
%    f = conformalmap(domain, range)
% will be able to select the correct method to create a conformal map,
% e.g.,
%    f = conformalmap(unitcircle, polygon).
% This second functionality has yet to be added.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

properties
  domain_                           % Region object.
  range_                            % Region object.
  function_list_
end

methods
  function f = conformalmap(domain, range, varargin)
    % Eventually, construct map given domain and range if called as base
    % class. Store domain and range otherwise.
    
    if ~nargin
      return
    end
    
    composition = false;
    if nargin >= 2 && all(cellfun(@(c) isa(c, 'conformalmap'), ...
                          [{domain, range}, varargin]))
      composition = true;
      argsok = true;
    elseif nargin == 2
      argsok = (isempty(domain) | isa(domain, 'region')) ...
               & (isempty(range) | isa(range, 'region'));
    else
      argsok = false;
    end
    if ~argsok
      error('CMT:InvalidArguments', ...
            'Expected domain and range region objects.')
    end
    
    if composition
      f.function_list_ = [{domain, range}, varargin];
      f.domain_ = f.function_list_{1}.domain;
      f.range_ =f.function_list_{end}.range;
    else
      f.domain_ = domain;
      f.range_ = range;
    end
  end
  
  function w = apply(f, z)
    % w = apply(f, z)
    %   Apply the conformal map to z.
    %
    % See the apply concept in the developer docs for more details.
    
    if iscomposition(f)
      % f is a composition, apply each map in turn.
      w = z;
      for k = 1:numel(f.function_list_)
        w = apply(f.function_list_{k}, w);
      end
      return
    end
    
    if ismethod(z, 'apply')
      % Try asking the target object first.
      try
        w = apply(z, f);
        return
      catch err
        if ~any(strcmp(err.identifier, ...
                       {'MATLAB:UndefinedFunction', 'CMT:NotDefined'}))
          rethrow(err)
        end
      end
    end
    
    % Try the apply defined by subclass.
    try
      if isa(z, 'gridcurves')
        w = cell(numel(z), 1);
        for k = 1:numel(w)
          w{k} = apply_map(f, z{k});
        end
        w = gridcurves(w);
      else
        w = apply_map(f, z);
      end
    catch err
      if strcmp(err.identifier, 'MATLAB:UndefinedFunction')
        msg = sprintf('Applying %s to %s is not defined.', class(f), class(z));
        if ~isempty(err.message)
          msg = sprintf('%s\n%s', msg, err.message);
        end
        error('CMT:NotDefined', msg)
      else
        rethrow(err)
      end
    end
  end
  
  function disp(f)
    fprintf('** conformalmap object **\n\n')
    if ~isempty(f.domain_)
      fprintf('---\nhas domain\n')
      disp(f.domain_)
    end
    if ~isempty(f.range_)
      fprintf('---\nhas range\n')
      disp(f.range_)
    end
    
    if iscomposition(f)
      fprintf('---\nthis map is a composition of:\n(in order of application)\n')
      for k = 1:numel(f.function_list_)
        disp(f.function_list_{k})
      end
    end
  end
  
  function d = domain(f)
    % Return map domain region.
    d = f.domain_;
  end
  
  function tf = iscomposition(f)
    tf = ~isempty(f.function_list_);
  end
  
  function f = mtimes(f1, f2)
    % Interpret f1*f2 as the composition f1(f2(z)).

    if ~(isa(f1, 'conformalmap') && isa(f2, 'conformalmap'))
      error('CMT:NotDefined', 'Operation is not defined.')
    end    
    f = conformalmap(f2, f1);
  end
  
  function out = plot(f)
    washold = ishold;
    
    if isempty(f.range_) || (~isempty(f.domain_) && isempty(grid(f.domain_)))
      warning('Map range not set or domain has an empty grid. No plot produced.')
      return
    end
    
    cah = newplot;
    hold on
    
    hg = plot(apply(f, grid(f.domain_)));
    hb = plot(f.range_);
    
    if ~washold
      plotdef.whitefigure(cah)
      axis(plotbox(f.range_))
      aspectequal
      axis off
      hold off
    end
    
    if nargout
      out = [hg; hb];
    end
  end
  
  function r = range(f)
    % Return map range region.
    r = f.range_;
  end
  
  function out = subsref(f, S)
    if length(S) == 1 && strcmp(S.type, '()')
      out = apply(f, S.subs{:});
    else
      out = builtin('subsref', f, S);
    end
  end
end

methods(Access=protected)
  function w = apply_map(~, z)
    % Identity map.
    w = z;
  end
end

end
