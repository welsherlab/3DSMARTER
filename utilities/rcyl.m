classdef rcyl < handle
    
    properties (SetAccess =immutable)
        P1
        P2
        Q
    end
    properties (SetAccess =protected)
        cyl % internal cyl obj
        
        theta_sign % apply sign when theta in cyl is obtained or used
        
        
        
        
    end
    methods
        function self = rcyl(p1,p2,q)
            self.cyl = cyl_obj(p1,p2,q);
            
            self.P1 = p1;
            self.P2 = p2;
            self.Q = q;
            
            if self.cyl.is_rh
                self.theta_sign = 1;
                
            else
                self.theta_sign = -1;
                
            end
            
            
        end
        
        
        function out = cart2cyl(self,cart_pos)
            % 210728: automatically add zeros in 3rd col for 2d inputs
            
            if length(cart_pos(1,:)) == 2
                cart_pos(:,3) = 0;
                %disp('added zeros in column 3 for cartesian coords')
            end
            
            out = self.cyl.cart2cyl_bat(cart_pos);
            out(:,2) = self.theta_sign.*out(:,2);
        end
        
        function out = cyl2cart(self,cyl_pos)
            
            if any(isnan(cyl_pos(:,2)))
                warning('nans in theta converted to 0s')
                cyl_pos( isnan(cyl_pos(:,2)),2) = 0;
            end
            
            cyl_pos(:,2) = self.theta_sign.*cyl_pos(:,2);
            out = self.cyl.cyl2cart_bat(cyl_pos);
        
        end
        
        function out = cart2sph(self,cart_pos)
            out = self.cyl2sph(self.cart2cyl(cart_pos));
            
        end
        function out = sph2cart(self,sph_pos)
            out = self.cyl2cart(self.sph2cyl(sph_pos));
            
        end        
        function out = cyl2sph(self,cyl_pos)
            r = sqrt(cyl_pos(:,1).^2 + cyl_pos(:,3).^2);
            theta = atan( cyl_pos(:,1)./cyl_pos(:,3))./pi.*180; %-90-90 degrees
            out = [abs(r) ,self.theta_sign.*wrapTo360( rcyl.theta_cyl2sph( cyl_pos(:,2) ) ) ,rcyl.atan_conv(theta)];
        end
        
        function out = sph2cyl(self,sph_pos)
            
            sph_pos(:,2) = wrapTo360(sph_pos(:,2));
            
            sph_pos(:,3) = abs(wrapTo180(sph_pos(:,3)));
            r = sph_pos(:,1).*sin(sph_pos(:,3)./180.*pi);
            theta= rcyl.theta_sph2cyl( sph_pos(:,2) );
            h = sph_pos(:,1).*cos(sph_pos(:,3)./180.*pi);
            
            out = [abs(r),self.theta_sign.*theta,h];
            
        end
        
    end
    
    methods (Static)
        function out = theta_sph2cyl(in) 
            %out = wrapTo180(in);
            %out = in;
            %out(in>180) = out(in>180) - 360;
            
            out = in - 180;
        end
        
        function out = atan_conv(in)
        out = in;
        out(in<0) = 180 +  in(in<0);
            
        end
        function out = theta_cyl2sph(in)
            %out = zeros(size(in));
            %out(in<0) = 360+in(in<0);
            %out(in>0) = in(in>0);
            
            % this is not true!
            %out = wrapTo360(in);
            
            out = in + 180;
        end
        
    end
    
end