classdef polargrid < zgrid
    
    properties
        radii = []
        angles = []
        type = 'mesh';
    end
    
    methods
        function g = polargrid(r,radii,angles,type)
            g = g@zgrid(r);
            if nargin==1
                radii = 8;
                angles = 12;
            end
            if length(radii)==1 && radii==round(radii)
                radii = (1:radii) / radii;
            end
            if length(angles)==1 && angles==round(angles)
                angles = 2*pi*(0:angles-1) / angles;
            end
            g.radii =  radii;
            g.angles = angles;
            
            if nargin > 3
                g.type = lower(type);
            end
            
            g = createsource(g);  % set up the initial grid
        end
        
        function g = apply(g,f)
            src = g.dataSource;
            switch(g.type)
                case 'curves'
                    for i = 1:length(g.radii)
                        img{1}{i} = apply(f,src{1}{i});
                    end
                    for i = 1:length(g.angles)
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

            switch(g.type)
                case 'curves'
                    for i = 1:length(src{1})
                        h1(i) = plot(img{1}{i}); hold on
                    end
                    for i = 1:length(src{2})
                        h2(i) = plot(img{2}{i});
                    end
                    colr = get(gca,'colororder');
                    set(h1,'color',colr(1,:))
                    set(h2,'color',colr(2,:))
                    h= [h1 h2];
                case 'mesh'
                    Z = img(:,[1:end 1]);
                    W = src(:,[1:end 1]);
                    h = PhasePlot.PhasePlot(Z,W,'e');
                otherwise
                    error('Unrecognized grid type.')
            end
            
            axis equal
            axis auto
            
            if ~washold
                hold off
            end
            
            if nargout > 0
                out = h;
            end
        end
    end
    
    methods (Hidden)
        function g = createsource(g)
            switch(g.type)
                case 'curves'
                    src = cell(1,2);
                    for i = 1:length(g.radii)
                        src{1}{i} = circle(0,g.radii(i));
                    end
                    for i = 1:length(g.angles)
                        src{2}{i} = segment(0,exp(1i*g.angles(i)));
                    end
                    g.dataSource = src;
                    g.dataImage = src;
                case 'mesh'
                    [R,T] = ndgrid(g.radii,g.angles);
                    g.dataSource = R.*exp(1i*T);
                    g.dataImage = g.dataSource;
                otherwise
                    error('Unrecognized grid type.')
            end
        end
         
    end
    
end
