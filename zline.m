classdef zline < curve
    
    properties
        point
    end
    
    methods
        function zl = zline(varargin)
            switch nargin
                case 1
                    zp = varargin{1};
                    % Look for two points in a vector.
                    if isa(zp,'double') && length(zp)==2 && all(~isinf(zp))
                        point_ = zp(1);
                        tangent_ = diff(zp);
                    else
                        error('Expected a vector of 2 finite points on the line.')
                    end
                case 2
                    % Look for a point and a tangent
                    [zp,ztan] = deal(varargin{:});
                    if isa(zp,'double') && isa(ztan,'double') && length(zp)==1 && length(ztan)==1
                        point_ = zp;
                        tangent_ = ztan;
                    else
                        error('Expected a point on the line and a tangent direction.')
                    end
                otherwise
                    error('Expected a vector of 2 points or a point and a tangent.')
            end
            
            zl = zl@curve(@position,@tangent,[-1 1]);
            
            function z = position(t)
                t = t./(1-t.^2);
                z = point_ + t*tangent_;
            end
            
            function zt = tangent(t)
                zt = tangent_;
            end
            
            zl.point = point_;
            
        end
        
        function str = char(zl)
            str = [ 'line passing through ',...
                num2str(zl.point),...
                ' and tangent to ',...
                num2str(zl.tangent(0)) ];
        end
        
        function d = dist(l,z)
            v = z - l.position(0);
            s = sign(l.tangent(0));
            d = abs(real(v)*real(s) + imag(v)*imag(s));
        end
        
        function tf = isinf(l)
            tf = true;
        end
        
        function l = truncate(l)
            l = segment(l.position(-0.5),l.position(0.5));
        end
    end
    
end   