classdef polargrid < zgrid
    
    % TODO: Allow translated center
    
    properties
        radii = []
        angles = []
        rbounds = [0 1]
        %tbounds = [-pi pi]    % TODO: Make this variable. (Need arcs.)
        type = 'mesh';
    end
    
    methods
        function g = polargrid(r,radii,angles,type,rbnd,tbnd)
            g = g@zgrid(r);
            
            if nargin > 3
                g.type = type;
                if nargin > 4
                    g.rbounds = rbnd;
%                     if nargin > 5
%                         g.tbounds = tbnd;
%                     end
                end
            end
            
            g.radii =  radii;
            g.angles = angles;
            g.PhasePlotType = 'e';
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
        
        function g = set.radii(g,r)
            rb = g.rbounds;                 % FIXME
            if length(r)==1 && r==round(r)
                if strcmp(g.type,'curves') && rb(1)==0  % FIXME
                    r = linspace(rb(1),rb(2),r+1);
                    r(1) = [];
                else
                    r = linspace(rb(1),rb(2),r);
                end
            end
            g.radii = r;
            g = createsource(g);
        end
        
        function g = set.angles(g,theta)
            if length(theta)==1 && theta==round(theta)
                theta = 2*pi*(0:theta-1) / theta;
            end
            g.angles = theta;
            g = createsource(g);
        end
        
        function g = set.type(g,type)
            g.type = type;
            g = createsource(g);
        end
        
%         function out = plot(g)
%             src = g.dataSource;
%             img = g.dataImage;
%             newplot
%             washold = ishold;
%             
%             pref = plotset;
%             plotargs = {'linewidth',pref.gridWidth,'color',pref.gridColor};
%             
%             switch(g.type)
%                 case 'curves'
%                     for i = 1:length(src{1})
%                         h1(i) = plot(img{1}{i}); hold on
%                     end
%                     for i = 1:length(src{2})
%                         h2(i) = plot(img{2}{i});
%                     end
%                     h = [h1 h2];
%                     set(h1,plotargs{:});
%                     set(h2,plotargs{:});
%                     
%                 case 'mesh'
%                     Z = img(:,[1:end 1]);
%                     W = src(:,[1:end 1]);
%                     h = PhasePlot.PhasePlot(Z,W,'e');
%                 otherwise
%                     error('Unrecognized grid type.')
%             end
%             
%             axis auto 
%             axis equal
%             
%             if ~washold
%                 hold off
%             end
%             
%             if nargout > 0
%                 out = h;
%             end
%         end
    end
    
    methods (Hidden)
        function g = createsource(g)
            switch(g.type)
                case 'curves'
                    src = cell(1,2);
                    for i = 1:length(g.radii)
                        src{1}{i} = circle(0,g.radii(i));
                    end
                    rb = g.rbounds;
                    for i = 1:length(g.angles)
                        tau = exp(1i*g.angles(i));                      
                        src{2}{i} = segment(rb(1)*tau,rb(2)*tau);
                    end
                    g.dataSource = src;
                    g.dataImage = src;
                case 'mesh'
                    [R,T] = ndgrid(g.radii,g.angles);
                    Z = R.*exp(1i*T);
                    if ~isempty(Z)
                        Z = Z(:,[1:end 1]);   % to join up in theta
                    end
                    g.dataSource = Z;          
                    g.dataImage = g.dataSource;
                otherwise
                    error('Unrecognized grid type.')
            end
        end
        
    end
    
end
