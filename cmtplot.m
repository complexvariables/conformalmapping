classdef cmtplot < cmtobject
% CMTPLOT collects common CMT plot tasks.
%
% This is probably more convoluted than it needs to be.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

properties(Constant)
    % colors
    black = 'k'
    grey = 0.75*[1, 1, 1]
    none = 'none'
end

methods
    function p = cmtplot
        % This is here to allow get/set functionality.

        get(p, plotset);
    end
end

methods(Static)
    function args = closedcurveargs()
        opts = get(cmtplot);
        args = { ...
            'color', opts.lineColor, ...
            'linewidth', opts.lineWidth, ...
            'linesmoothing', opts.lineSmoothing ...
        };
    end

    function args = fillargs()
        args = {cmtplot.fillcolor, 'edgecolor', cmtplot.filledgecolor};
    end

    function c = fillcolor()
        c = cmtplot.grey;
    end

    function c = filledgecolor()
        c = cmtplot.none;
    end

    function args = gridargs()
        opts = get(cmtplot);
        args = {
            'color', opts.gridColor, ...
            'linewidth', opts.lineWidth, ...
            'linesmoothing', opts.lineSmoothing ...
        };
    end

    function [args, gargs] = pullgridargs(arglist)
        % Separate grid args from cell array arglist. Must have an even 
        % number of entries of the {'name', value} pair form.

        n = numel(arglist);
        if mod(n, 2)
            error('CMT:InvalidArgument', 'Expected name/value pairs.')
        end

        recognized = properties(gridset);
        gargs = {};
        idx = 1:n;
        for k = 1:2:n-1
            match = strcmpi(arglist{k}, recognized);
            if ~any(match)
                continue
            end
            gargs(numel(gargs)+(1:2)) = ...
                {recognized(match), arglist{k+1}}; %#ok<AGROW>
            idx = idx(idx ~= k & idx ~= k+1);
        end
        if isempty(idx)
            args = {};
        else
            args = arglist(idx);
        end
    end

    function whitefigure(fig)
        if ~nargin
            fig = gcf;
        end
        set(fig, 'color', 'white');
    end
end
    
end
