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
    if ~nargin
      return
    end
    
    argsok = false;
    composition = false;
    anonymous = false;
    if nargin >= 2
      if all(cellfun(@(c) isa(c, 'conformalmap'), ...
                          [{domain, range}, varargin]))
        composition = true;
        argsok = true;
      else
        argsok = (isempty(domain) | isa(domain, 'region')) ...
                 & (isempty(range) | isa(range, 'region'));
        if nargin == 3
          anonymous = isa(varargin{1}, 'function_handle');
          argsok = argsok & anonymous;
        elseif nargin > 3
          argsok = false;
        end
      end
    end
    if ~argsok
      error('CMT:InvalidArgument', ...
            'Expected domain and range region objects.')
    end
    
    if composition
      f.function_list_ = [{domain, range}, varargin];
      f.domain_ = f.function_list_{1}.domain;
      f.range_ = f.function_list_{end}.range;
    else
      f.domain_ = domain;
      f.range_ = range;
      if anonymous
        f.function_list_ = varargin;
      end
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
  
  function tf = isanonymous(f)
    tf = numel(f.function_list_) == 1 ...
         && isa(f.function_list_{1}, 'function_handle');
  end
  
  function tf = iscomposition(f)
    tf = numel(f.function_list_) > 1;
  end
  
  function f = mtimes(f1, f2)
    % Interpret f1*f2 as the composition f1(f2(z)).

    if ~(isa(f1, 'conformalmap') && isa(f2, 'conformalmap'))
      error('CMT:NotDefined', 'Operation is not defined.')
    end    
    f = conformalmap(f2, f1);
  end
  
  function out = plot(f, varargin)
    washold = ishold;
    
    if isempty(f.range_) || (~isempty(f.domain_) && isempty(grid(f.domain_)))
      warning('Map range not set or domain has an empty grid. No plot produced.')
      return
    end
    
    cah = newplot;
    hold on

    [pargs, gargs] = plotdef.pullgridargs(varargin);
    hg = plot(apply(f, grid(f.domain_, gargs{:})));
    hb = plot(f.range_, pargs{:});
    
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
  
  function varargout = subsref(f, S)
    if numel(S) == 1 && strcmp(S.type, '()')
      [varargout{1:nargout}] = apply(f, S.subs{:});
    else
      [varargout{1:nargout}] = builtin('subsref', f, S);
    end
  end
end

methods(Access=protected)
  function w = apply_map(f, z)
    if isanonymous(f)
      w = f.function_list_{1}(z);
    else
      % Default map is identity.
      w = z;
    end
  end
end

end
