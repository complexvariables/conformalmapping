classdef conformalmap < cmtobject
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
    theDomain                           % Region object.
    theRange                            % Region object.
    functionList
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
            f.functionList = [{domain, range}, varargin];
            f.theDomain = f.functionList{1}.domain;
            f.theRange = f.functionList{end}.range;
        else
            f.theDomain = domain;
            f.theRange = range;
            if anonymous
                f.functionList = varargin;
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
            for k = 1:numel(f.functionList)
                w = apply(f.functionList{k}, w);
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
                    w{k} = applyMap(f, z{k});
                end
                w = gridcurves(w);
            else
                w = applyMap(f, z);
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
        if ~isempty(f.theDomain)
            fprintf('---\nhas domain\n')
            disp(f.theDomain)
        end
        if ~isempty(f.theRange)
            fprintf('---\nhas range\n')
            disp(f.theRange)
        end

        if iscomposition(f)
            fprintf('---\nthis map is a composition of:\n(in order of application)\n')
            for k = 1:numel(f.functionList)
                disp(f.functionList{k})
            end
        end
    end

    function d = domain(f)
        % Return map domain region.
        d = f.theDomain;
    end

  function w = evaluate(f, z)
    % w = apply(f, z)
    %   Apply the conformal map to z.
    %
    % See the apply concept in the developer docs for more details.
    
    if iscomposition(f)
      % f is a composition, apply each map in turn.
      w = z;
      for k = 1:numel(f.function_list_)
        w = evaluate(f.function_list_{k}, w);
      end
    else   
        % Try asking the target object first.
        try
            w = evaluate(z, f);
            return
        catch err
            if ~any(strcmp(err.identifier, ...
                    {'MATLAB:UndefinedFunction', 'CMT:NotDefined'}))
                rethrow(err)
            end
        end
    end
 
  end

    function tf = isanonymous(f)
        tf = numel(f.functionList) == 1 ...
            && isa(f.functionList{1}, 'function_handle');
    end

    function tf = iscomposition(f)
        tf = numel(f.functionList) > 1;
    end

    function out = plot(f, varargin)
        %FIXME: Some maps are going to want the grid to be in the range,
        %not the domain.
        washold = ishold;

        if ~isempty(f.theDomain) && ~hasgrid(f.theDomain)
            error('Domain has an empty grid. No plot produced.')
        end

        cah = newplot;
        hold on
        
        % Separate grid construction and plot 'name'/value pairs.
        [gargs, pargs] = separateArgs(get(f.theDomain), varargin{:}); 
        hg = plot(apply(f, grid(f.theDomain, gargs{:})), pargs{:});
        if ~isempty(f.theRange)
            hb = plot(boundary(f.theRange), pargs{:});
            pb = plotbox(f.theRange);
        else
            pb = axis;
        end

        if ~washold
            cmtplot.whitefigure(get(cah, 'parent'))
            axis(pb)
            aspectequal
            if cmtplot.hasNewGraphics
                % Bug in new graphics won't show smoothed lines if this
                % isn't done when drawing on a figure created by the plot
                % call. Just need some way to force a redraw after
                % the lines are put up.
                drawnow
            end
            axis off
            hold off
        end

        if nargout
            out = [hg; hb];
        end
    end

    function r = range(f)
        % Return map range region.
        r = f.theRange;
    end

    function varargout = subsref(f, S)
        if numel(S) == 1 && strcmp(S.type, '()')
            [varargout{1:nargout}] = apply(f, S.subs{:});
        else
            [varargout{1:nargout}] = builtin('subsref', f, S);
        end
    end
end

methods
% Arithmetic operators.
    function f = mtimes(f1, f2)
        % Interpret f1*f2 as the composition f1(f2(z)), or as scalar times
        % the image.

        % Make sure a conformalmap is first.
        if isnumeric(f1)
            tmp = f1;  f1 = f2;  f2 = tmp;
        end

        if isnumeric(f2)
            const = f2;  % for naming
            f = conformalmap(f1,@(z) const*z);
        elseif isa(f2,'conformalmap')
            % Do a composition of the two maps.
            f = conformalmap(f2, f1);
        else
            error('CMT:conformalmap:mtimes',...
                'Must either multiply by a scalar or compose two maps.')
        end
    end

    function g = minus(f,c)
        if isnumeric(f)
            g = plus(-f,c);
        else
            g = plus(f,-c);
        end
    end

    function g = plus(f,c)
        % Defines adding a constant (translation of the range).

        if isnumeric(f)
            tmp = c;  c = f;  f = tmp;
        end

        if ~isnumeric(c) || (length(c) > 1)
            error('CMT:conformalmap:plus',...
                'You may only add a scalar to a conformalmap.')
        end

        g = conformalmap( f, @(z) z+c );
    end

    function g = mpower(f,n)

        if (n ~= round(n)) || (n < 0)
            error('CMT:conformalmap:integerpower',...
                'Power must be a positive integer.')
        end

        g = f;
        p = floor(log2(n));
        for k = 1:p
            g = conformalmap(g,g);
        end
        for k = n-p
            g = conformalmap(g,f);
        end
    end

    function g = uminus(f)
        g = conformalmap( f, @(z) -z );
    end
end

methods(Access=protected)
    function w = applyMap(f, z)
        if isanonymous(f)
            w = f.functionList{1}(z);
        else
            % Default map is identity.
            if ~isa(f, 'conformalmap')
                warning('CMT:BadThings', ...
                    ['Identity map used when conformalmap is subclassed.\n' ...
                    '(Define applyMap in subclass %s?)'], class(f))
            end
            w = z;
        end
    end
end

end
