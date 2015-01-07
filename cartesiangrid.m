classdef cartesiangrid < zgrid
    
    properties
        real = []
        imag = []
        bounds = [-1 1 -1 1];
        type = 'mesh';
    end
    
    methods
        function g = cartesiangrid(r,bounds_,real_,imag_,type_)
            g = g@zgrid(r);
            
            if nargin > 4
                g.type = type_;
            end
            
            g.bounds = bounds_;
            g.real = real_;
            g.imag = imag_;
            
            g = createsource(g);  % set up the initial grid

        end
        
        function g = apply(g,f)
            src = g.dataSource;

            switch(g.type)
                case 'curves'
                    for i = 1:length(src{1})
                        img{1}{i} = apply(f,src{1}{i});
                    end
                    for i = 1:length(src{2})
                        img{2}{i} = apply(f,src{2}{i});
                    end
                    g.dataImage = img;
                case 'mesh'
                    g.dataImage = apply(f,g.dataSource);
                otherwise
                    error('Unrecognized grid type.')
            end
        end
        
        function out = plot(g)
            src = g.dataSource;
            img = g.dataImage;
            
            newplot
            washold = ishold;
            hold on
            
            pref = cmtgetpref('grid');
            plotargs.one = {'linewidth',pref.curvewidth,'color',pref.curvecolor(1,:)};
            plotargs.two = {'linewidth',pref.curvewidth,'color',pref.curvecolor(2,:)};

            switch(g.type)
                case 'curves'
                    for i = 1:length(src{1})
                        h1(i) = plot(img{1}{i});
                    end
                    for i = 1:length(src{2})
                        h2(i) = plot(img{2}{i});
                    end
                    h = [h1 h2];
                    set(h1,plotargs.one{:});
                    set(h2,plotargs.two{:});
                case 'mesh'
                    Z = img(:,[1:end 1]);
                    W = src(:,[1:end 1]);
                    h = PhasePlot.PhasePlot(Z,W,'v');
                otherwise
                    error('Unrecognized grid type.')
            end
            
            axis auto
            axis equal
            if ~washold
                hold off
            end
            
            if nargout > 0
                out = h;
            end
           
        end
    
    function g = set.imag(g,y)
        if length(y)==1 && y==round(y)
            y = equidistribute(g,y,3:4);
        end
        g.imag = y;
        g = createsource(g);
    end
    
    function g = set.real(g,x)
        if length(x)==1 && x==round(x)
            x = equidistribute(g,x,1:2);
        end
        g.real = x;
        g = createsource(g);
    end

    
    end
    methods (Hidden)
        function g = createsource(g)
            % Set up the source grid with whichever type of implementation
            % is selected. 
            b = g.bounds;
            switch(g.type)
                case 'curves'
                    src = cell(1,2);
                    % For vertical segments with infinite endpoints, have
                    % to get the tangent direction right. 
                    for i = 1:length(g.real)
                        x = g.real(i);
                        if isinf(b(3))
                            left = homog(-1i,0);
                        else
                            left = x+1i*b(3);
                        end
                        if isinf(b(4))
                            right = homog(1i,0);
                        else
                            right = x+1i*b(4);
                        end
                        src{1}{i} = segment(left,right);
                    end
                    for i = 1:length(g.imag)
                        y = g.imag(i);
                        src{2}{i} = segment(b(1)+1i*y,b(2)+1i*y);
                    end
                    g.dataSource = src;
                    g.dataImage = src;
                case 'mesh'
                    [X,Y] = ndgrid(g.real,g.imag);
                    g.dataSource = complex(X,Y);
                    g.dataImage = g.dataSource;
                otherwise
                    error('Unrecognized grid type.')
            end
        end
    end
    
    methods (Access=private)
        % Used to "equally" distribute values in a
        % finite/(semi)infinite interval.
        function q = equidistribute(g,n,indx)
            bounds = g.bounds(indx);
            if isinf(bounds(1)) && isinf(bounds(2))
                t = linspace(-1,1,n+2);
                t = t(2:n+1);
                q = 2*t./(1-t.^2);
            elseif isinf(bounds(1))
                t = linspace(-1,0,n+1);
                t = t(2:n+1);
                q = t./(1+t) + bounds(2);
            elseif isinf(bounds(2))
                t = linspace(0,1,n+1);
                t = t(1:n);
                q = t./(1-t) + bounds(1);
            else
                q = linspace(bounds(1),bounds(2),n);
            end
        end
    end
    
end
