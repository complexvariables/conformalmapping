classdef segment < curve

    properties
        endpoints  
    end
    
    methods
        function c = segment(za,zb)
            if nargin==0
                % Return an abstract segment with unknown geometry.
                return
            end
            
            if isinf(za) && isinf(zb)
                error('At least one endpoint must be finite.')
            elseif isinf(zb)
                delta = exp(1i*angle(zb));
                position = @(t) za + t./(1-t) * delta;
            elseif isinf(za)
                delta = -exp(1i*angle(za));
                position = @(t) zb - (1-t)./t * delta;
            else
                delta = zb-za;
                position = @(t) za + t * delta;
            end
            
            c = c@curve(position,@(t) delta*ones(size(t)),[0 1]);
            c.endpoints = [za zb];

        end
        
        function disp(s)
            fprintf('line segment with endpoints: ')
            disp(s.endpoints)
        end
        
        function d = dist(s, z)
            % Distance between point and line segment.
            error('TODO: Not supported yet.')
        end
        
        function z = intersect(s1,s2)
            % Calculate segment intersection.
            error('TODO: Not supported yet.')            
       end
        
        function tf = isinf(s)
            tf = any(isinf(s.endpoints));
        end
        
         function c = uminus(c)
            c = uminus@curve(c);
            c.endpoints = -c.endpoints;
        end
        
        function c = plus(c,z)
            c = plus@curve(c,z);
            c.endpoints = c.endpoints + z;
        end

        function c = mtimes(c,z)
            c = mtimes@curve(c,z);
            c.endpoints = c.endpoints * z;
        end
     
   end
    
    
end
