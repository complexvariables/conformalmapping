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
    ver84 = [ 0   4.4700e-01   7.4100e-01;...
        8.5000e-01   3.2500e-01   9.8000e-02;...
        9.2900e-01   6.9400e-01   1.2500e-01;...
        4.9400e-01   1.8400e-01   5.5600e-01;...
        4.6600e-01   6.7400e-01   1.8800e-01;...
        3.0100e-01   7.4500e-01   9.3300e-01;...
        6.3500e-01   7.8000e-02   1.8400e-01];
end

methods
    function p = cmtplot
        % This is here to allow get/set functionality.

        get(p, plotset);
    end
end

methods(Static)
    function [args, exargs] = closedcurveArgs(varargin)
        opts = get(cmtplot);
        if nargout > 1
            [opts, exargs] = set(opts, varargin{:});
        else
            opts = set(opts, varargin{:});
        end
        args = { ...
            'color', opts.lineColor, ...
            'linewidth', opts.lineWidth ...
        };
        if ~cmtplot.hasNewGraphics
            args(5:6) = {'linesmoothing', opts.lineSmoothing};
        end
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

    function [args, exargs] = gridArgs(varargin)
        opts = get(cmtplot);
        if nargout < 2
            opts = set(opts, varargin{:});
        else
            [opts, exargs] = set(opts, varargin{:});
        end
        args = {
            'color', opts.gridColor, ...
            'linewidth', opts.lineWidth ...
        };
        if ~cmtplot.hasNewGraphics
            args(5:6) = {'linesmoothing', opts.lineSmoothing};
        end
    end
    
    function tf = hasNewGraphics()
        tf = ~verLessThan('matlab', '8.4');
    end
    
    function tf = isFigHandle(val)
        if exist('isgraphics', 'builtin') == 5
            tf = isgraphics(val, 'figure');
        else
            % Revert to older way.
            if ishghandle(tmp)
            	tf = strcmp(get(tmp, 'type'), 'figure');
            else
                tf = false;
            end
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
