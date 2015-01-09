classdef carlesongrid < zgrid
    
    properties
        levels
    end
       

    methods
        function g = carlesongrid(R)
            %TODO: Accept options.
            g = g@zgrid(R);
            g.levels = get(gridset,'numLevels');
            
            nu = 32; % Base radial line number.
            r = 0.6*R.radius; % Base circle radius.

            gc = cell(2 + g.levels + 2^(g.levels-1)*nu, 1);

            % Level 0 circle.
            gc{1} = circle(R.center,r); 

            idx = 1;
            for j = 1:g.levels
                if j > 1
                    nuj = 2^(j-2)*nu;
                else
                    nuj = nu;
                end
                dt = 2*pi/nuj;
                offset = (j > 1)*dt/2;
                for k = 1:nuj
                    endpts = R.center + ...
                        R.radius*[r 1-1e-8]*exp(1i*(offset + (k-1)*dt));
                    gc{idx + k} = segment(endpts(1),endpts(2));
                end

                idx = idx + nuj + 1;
                r = (1 + r)/2;
                gc{idx} = circle(R.center,r*R.radius); 
            end
            gc{end} = circle(R.center,R.radius);
            
            g.dataSource = gc;
            g.dataImage = gc;
        end
        
        function g = apply(g,f)
            src = g.dataSource;
            for i = 1:length(src)
                img{i} = apply(f,src{i});
            end
            g.dataImage = img;
        end
                
        function out = plot(g)
            src = g.dataSource;
            img = g.dataImage;
            newplot
            washold = ishold;
            
            pref = plotset;
            plotargs = {'linewidth',pref.gridWidth,'color',pref.gridColor};
            
            for i = 1:length(src)
                h(i) = plot(img{i}); hold on
            end
            set(h,plotargs{:});
            
            if ~washold
                axis auto
                axis equal
                hold off
            end
            
            if nargout > 0
                out = h;
            end
        end
    end
        
end
