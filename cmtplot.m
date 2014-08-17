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
    function args = closedcurveargs(varargin)
        opts = get(cmtplot);
        opts = set(opts, varargin{:});
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

    function args = gridargs(varargin)
        opts = get(cmtplot);
        opts = set(opts, varargin{:});
        args = {
            'color', opts.gridColor, ...
            'linewidth', opts.lineWidth, ...
            'linesmoothing', opts.lineSmoothing ...
        };
    end

    function [gargs, args] = pullGridArgs(varargin)
        % Separate grid args from cell array arglist. Must have an even 
        % number of entries of the {'name', value} pair form.

        if mod(nargin, 2)
            error('CMT:InvalidArgument', 'Expected name/value pairs.')
        end

        gargs = {};
        idx = 1:nargin;
        for k = 1:2:nargin-1
            match = regexpi(varargin{k}, '^grid.*$', 'match');
            if numel(match) == 0
                continue
            end
            match = match{1};
            gargs(numel(gargs)+(1:2)) = ...
                {match, varargin{k+1}}; %#ok<AGROW>
            idx = idx(idx ~= k & idx ~= k+1);
        end
        if isempty(idx)
            args = {};
        else
            args = varargin(idx);
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
