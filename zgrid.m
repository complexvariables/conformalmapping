classdef (Abstract) zgrid
    % GRIDCURVES class holds grid curves for regions.
    %
    % grid = gridcurves(curves)
    % The given cell array curves is stored in the object. This is mainly to
    % facilitate the use of plot(grid) to give a standard look to grid plots.
    %
    % Stores grid curves as entries in cell array. Consider
    %    gd = gridcurves;
    % We overload gd(n) and gd(n,m) to retrieve those cell array entries. Why
    % not just use the '{}' syntax? Wouldn't it be clearer we're using cell arrays?
    
    % This file is a part of the CMToolbox.
    % It is licensed under the BSD 3-clause license.
    % (See LICENSE.)
    
    % Copyright Toby Driscoll, 2014.
    % (Re)written by Everett Kropf, 2014,
    % adapted from an idea by Toby Driscoll, 20??.
    
    properties
        region
        PhasePlotType = 'e';
    end
    
    properties (SetAccess = protected)
        dataSource = {}
        dataImage = {}
    end
    
    methods
        function g = zgrid(r)
            if ~nargin
                return
            end
            
            g.region = r;
        end
        
        function s = char(g)
            s = sprintf('grid object (%s) in:\n',class(g));
            s = [s, '   ',char(g.region)];
            s = strrep(s,'\n','\n   ');
            s = deblank(s);
        end
        
        function disp(g)
            fprintf(char(g))
        end
        
        function out = plot(g)
            src = g.dataSource;
            img = g.dataImage;
            
            %newplot
            washold = ishold;
            hold on
            
            pref = plotset;
            plotargs = {'linewidth',pref.gridWidth,'color',pref.gridColor};
            
            % TODO: Base plot type on the data type rather than what the
            % grid says it is? Might simplify life with Carleson.
            switch(g.type)
                case 'curves'
                    h1 = [];
                    for i = 1:length(src{1})
                        h1(i) = plot(img{1}{i});
                    end
                    h2 = [];
                    for i = 1:length(src{2})
                        h2(i) = plot(img{2}{i});
                    end
                    h = [h1 h2];
                    set(h1,plotargs{:});
                    set(h2,plotargs{:});
                case 'mesh'
                    Z = img;
                    W = src;
                    y = g.PhasePlotType;
                    h = PhasePlot.PhasePlot(Z,W,y);
                otherwise
                    error('Unrecognized grid type.')
            end
            
            if ~washold
                axis auto
                axis equal
                box on
                hold off
            end
            
            if nargout > 0
                out = h;
            end
            
        end
        
        function out = rsplot(g)
            src = g.dataSource;
            img = g.dataImage;
            
            %newplot
            washold = ishold;
            
            pref = plotset;
            plotargs = {'linewidth',pref.gridWidth,'color',pref.gridColor};
            
            % TODO: Base plot type on the data type rather than what the
            % grid says it is? Might simplify life with Carleson.
            switch(g.type)
                case 'curves'
                    for i = 1:length(src{1})
                        h1(i) = rsplot(img{1}{i});
                        hold on
                    end
                    for i = 1:length(src{2})
                        h2(i) = rsplot(img{2}{i});
                    end
                    h = [h1 h2];
                    set(h1,plotargs{:});
                    set(h2,plotargs{:});
                otherwise
                    error('Unrecognized grid type.')
            end
            
            if ~washold
                axis equal
                hold off
            end
            
            if nargout > 0
                out = h;
            end
            
        end
        
    end
    
    methods (Abstract)
        w = apply(g,f)    % apply function f to grid g
    end
    
end


