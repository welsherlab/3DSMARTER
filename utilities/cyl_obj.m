classdef cyl_obj <handle
    
    properties(SetAccess = immutable)
        theta_offset % this is important to convert cartesian coords to cyl coords
        Q_1 % begin point of center axis
        Q_2 % end point of center axis, this helps decide the direction of z' (h) axis
        ref_pt % ref point that defines theta = 0 
        params % thee params (A,B,C) that describes a plane determined by Ax+By+Cz+1 = 0 or Ax+By+Cz = 0
        d_0 % distance between Q1 and Q2
        issing % is the plane determined by Q1, Q2 and point0 singular? (is [0,0,0] on the plane?)
        orig_angle_sign % original angle sign
        
        % added 210703
        k1 
        b1
        k2
        b2
        k3
        b3
    end
    
    properties
        angle_sign
    end
    
    methods
        function  self = cyl_obj(Q1,Q2,point0)
            
            if ~isvector(Q1)||~isvector(Q2)||~isvector(point0)||length(Q1)~=3 ... 
                    ||length(Q2)~=3||length(point0)~=3
                error('inputs must be vectors and have 3 elements!')
            end
            
            Q1 = Q1(:)'; Q2 = Q2(:)'; point0 = point0(:)';
            
            if Q1 == Q2
               error('Q1 and Q2 must be different points!') 
            end
            
            if cross(Q1-point0,Q2-point0)==0
                error('refence point cannot be on the line defined by Q1 and Q2!')
            end
            
            A = zeros(4);
            A(:,4) = 1;
            A(1,1:3) = Q1;
            A(2,1:3) = Q2;
            A(3,1:3) = point0;
  
            self.Q_1 = Q1;
            self.Q_2 = Q2;
            self.ref_pt = point0;
            self.d_0 = sqrt(sum((Q2(:) - Q1(:)).^2));
            
            normal = Q1 - Q2;
            v=null(normal);
            center = Q1;
            radius = sqrt(sum((Q1 - Q2).^2));
            
            f = @(theta) sum((repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta)) - point0').^2);
            self.theta_offset = fminsearch(f,0);
            
            left_mat = [Q1;Q2;point0];
            
            if det(A) == 0
               self.issing = 1;
               %disp('note: [0,0,0] is on the plane determined by the given 3 points')
               self.params = cross(Q1 - point0, Q2 - point0);
            else
                self.issing = 0;
                self.params = linsolve(left_mat,-ones(3,1)); 
            end
            [self.k1,self.b1,~] = lf_obj.linear_fit_single([0,self.d_0],[self.Q_1(1),self.Q_2(1)]);
            [self.k2,self.b2,~] = lf_obj.linear_fit_single([0,self.d_0],[self.Q_1(2),self.Q_2(2)]);
            [self.k3,self.b3,~] = lf_obj.linear_fit_single([0,self.d_0],[self.Q_1(3),self.Q_2(3)]);            
            self.angle_sign = 1;
            self.orig_angle_sign = 1;
            test_pts = [1.111,2.345,3.555; 19.94,4,24; 13,82,86];
            for i = 1:length(test_pts(:,1))
                curr_pts = test_pts(i,:);
                curr_pts_post_conv = self.cyl2cart(self.cart2cyl(curr_pts));
                if get_dist(curr_pts,curr_pts_post_conv) > 1
%                     disp(['input of ' num2str(curr_pts) ' is converted to ' num2str(curr_pts_post_conv) ... 
%                         ' in a test run, which exceeds threshold, and thus '  'angle_sign is set to -1'])
                    self.orig_angle_sign = -1;
                    self.angle_sign = -1;
                    break,
                end
                
            end
            

             
            
        end
       
        function cyl_coords = cart2cyl(self,test_pt) % returns cyl coords (R,theta(in degrees -180 to +180),z' (or h))
            Q1 = self.Q_1; Q2 = self.Q_2; point0 = self.ref_pt; pl_params = self.params;
            
            point0_proj = pt_proj_on_line(Q1,Q2,point0);
            test_pt_proj = pt_proj_on_line(Q1,Q2,test_pt);
            vec1 = point0(:) - point0_proj(:);
            vec2 = test_pt(:) - test_pt_proj(:);
            
            if self.issing == 0
            test_pt_sign =  sign(  pl_params(1).*(test_pt(1)) + pl_params(2).*(test_pt(2)) + pl_params(3).*(test_pt(3)) + 1); 
            else
                %warning('matrix is singular')
            test_pt_sign = sign(  pl_params(1).*(test_pt(1)) + pl_params(2).*(test_pt(2)) + pl_params(3).*(test_pt(3)) + 0);    
            end
            test_pt_rad = test_pt_sign*acos(dot(vec1,vec2)./( norm(vec1).*norm(vec2) )); 
            test_pt_deg = test_pt_rad./pi.*180;
            test_pt_r = sqrt(sum((test_pt(:) - test_pt_proj(:)).^2));
            
            d0 = sqrt(sum((Q2(:) - Q1(:)).^2));
            d1 = sqrt(sum((test_pt_proj(:) - Q1(:)).^2));
            d2 = sqrt(sum((test_pt_proj(:) - Q2(:)).^2));
            test_pt_z = -d1.*sign([d2>d1 && d2>d0]-.5);
            
            % added 210724\
            if test_pt_sign == 0 && acos(dot(vec1,vec2)./( norm(vec1).*norm(vec2) ))>3.13
                test_pt_deg = 180;
                
            end
            
            cyl_coords = [test_pt_r,test_pt_deg,test_pt_z]; % R,theta,z
            %error('Even fortune cookies are getting cynical.')
        end
        
        function cyl_coords = cart2cyl_bat(self,test_pt) % added 210423
            
            Q1 = self.Q_1; Q2 = self.Q_2; point0 = self.ref_pt; pl_params = self.params;
            point0_proj = pt_proj_on_line(Q1,Q2,point0)';
            
            len = cyl_obj.check_nby3(test_pt);
            cyl_coords = zeros(len,3);
            
            P = repmat(self.Q_1,[len,1]);
            Q = repmat(self.Q_2,[len,1]);
            
            % get r first, seems simple
            projs= cyl_obj.pt2line(P,Q,test_pt); % project points
            r = sqrt(sum((test_pt - projs).^2,2));
            
            % next, get z'
            
            zp0 = sqrt(sum((P - projs).^2,2)); % this may not have correct signs!
            zp1 = sqrt(sum((Q - projs).^2,2)); % this is the distance to Q
            
            % this will return an error, && must return a scalar!
            %tf = (zp1>self.d_0) && (abs(zp1-zp0-self.d_0)<1e-9*self.d_0);
            
           % 210724£ºfollowing line removed, I don't know what this is doing 
            %tf = and(zp1>self.d_0,abs(zp1-zp0-self.d_0)<1e-9*self.d_0);
            % the right way to do this is, sign is -1 only when zp2-zp1 =
            % d0
            
            % added this line 210724
            tf = abs(zp1-zp0-self.d_0)<1e-5*self.d_0 ;
            
            zp = -sign(double(tf)-.5).*zp0;
            
            clearvars zp0 zp1 tf
            
            % get theta
            if self.issing == 0
            test_pt_sign = sign(  pl_params(1).*(test_pt(:,1)) + pl_params(2).*(test_pt(:,2)) + pl_params(3).*(test_pt(:,3)) + 1);
            else
            test_pt_sign = sign(  pl_params(1).*(test_pt(:,1)) + pl_params(2).*(test_pt(:,2)) + pl_params(3).*(test_pt(:,3)) + 0);    
            end           
            
            vec1 = repmat(point0-point0_proj,[len,1]);
            vec2 = test_pt - projs;
            
            clearvars projs
            
            
            norm_v1 = sqrt( sum(vec1.^2,2 ));
            norm_v2 = sqrt( sum(vec2.^2,2 ));
            test_rad = test_pt_sign.*acos(dot(vec1,vec2,2)./(norm_v1.*norm_v2));
            
           tmp_idx = [];
            % added 210724
            if any(test_pt_sign == 0)
            
                
                disp([num2str(length(find(test_pt_sign==0))) ' of points on the ref plane found' ])
                
                tmp_idx = intersect( find(dot(vec1,vec2,2)./(norm_v1.*norm_v2)<-.99),find(test_pt_sign==0) );
                if ~isempty(tmp_idx)
                    warning([num2str(length(tmp_idx)) ' of points with degree = 180' ])
                end
                
            end
   
            
            clearvars vec1 vec2 norm_v1 norm_v2
            final_deg= test_rad./pi.*180;
            final_deg(tmp_idx) = 180; % added 210724
            
            cyl_coords(:,1) = r;
            cyl_coords(:,2) = final_deg;
            cyl_coords(:,3) = zp;
            

            
        end
        
        function cyl_coords = cart2cyl_bat_gpu(self,test_pt) % added 210423
            
            Q1 = self.Q_1; Q2 = self.Q_2; point0 = self.ref_pt; pl_params = self.params;
            point0_proj = pt_proj_on_line(Q1,Q2,point0)';
            
            len = cyl_obj.check_nby3(test_pt);
            cyl_coords = gpuArray(zeros(len,3));
            
            P = repmat(self.Q_1,[len,1]);
            Q = repmat(self.Q_2,[len,1]);
            
            % get r first, seems simple
            projs= gpuArray(cyl_obj.pt2line(P,Q,test_pt)); % project points
            r = gpuArray(sqrt(sum((test_pt - projs).^2,2)));
            
            % next, get z'
            
            zp0 = gpuArray(sqrt(sum((P - projs).^2,2))); % this may not have correct signs!
            zp1 = gpuArray(sqrt(sum((Q - projs).^2,2))); % this is the distance to Q
            
            % this will return an error, && must return a scalar!
            %tf = (zp1>self.d_0) && (abs(zp1-zp0-self.d_0)<1e-9*self.d_0);
            
            tf = and(zp1>self.d_0,abs(zp1-zp0-self.d_0)<1e-9*self.d_0);
            zp = sign(double(tf)+.5).*zp0;
            
            clearvars zp0 zp1 tf
            
            % get theta
            if self.issing == 0
            test_pt_sign = sign(  pl_params(1).*(test_pt(:,1)) + pl_params(2).*(test_pt(:,2)) + pl_params(3).*(test_pt(:,3)) + 1);
            else
            test_pt_sign = sign(  pl_params(1).*(test_pt(:,1)) + pl_params(2).*(test_pt(:,2)) + pl_params(3).*(test_pt(:,3)) + 0);    
            end           
            
            vec1 = gpuArray(repmat(point0-point0_proj,[len,1]));
            vec2 = gpuArray(test_pt - projs);
            
            clearvars projs
            
            
            norm_v1 = sqrt( sum(vec1.^2,2 ));
            norm_v2 = sqrt( sum(vec2.^2,2 ));
            test_rad = test_pt_sign.*acos(dot(vec1,vec2,2)./(norm_v1.*norm_v2));
            
            clearvars vec1 vec2 norm_v1 norm_v2
            final_deg= test_rad./pi.*180;
            
            cyl_coords(:,1) = r;
            cyl_coords(:,2) = final_deg;
            cyl_coords(:,3) = zp;
            
            cyl_coords = gather(cyl_coords);
            
            %error('!!')
            
        end
        
        
        function cart_coords = cyl2cart(self,cyl_coords) % return cartesian coords
            curr_r = cyl_coords(1);
            curr_deg = self.angle_sign*cyl_coords(2);
            curr_z = cyl_coords(3);
            d0 = self.d_0;
            Q1 = self.Q_1;
            Q2 = self.Q_2;
            
            %[k1,b1,~] = lf_obj.linear_fit_single([0,d0],[Q1(1),Q2(1)]);
            %[k2,b2,~] = lf_obj.linear_fit_single([0,d0],[Q1(2),Q2(2)]);
            %[k3,b3,~] = lf_obj.linear_fit_single([0,d0],[Q1(3),Q2(3)]);
            
            center = [self.k1*curr_z+self.b1,self.k2*curr_z+self.b2,self.k3*curr_z+self.b3];
            
            % I don't know what these commented codes are used for but they
            % clearly don't work 210528
            %d1 = abs(curr_z);
            %d2 = abs(d1 -sign(curr_z)*d0);
            %center = (d1./d2.*Q2 - sign((d1>d0||d2>d0)-.5)* Q1)./(d1/d2 - sign((d1>d0||d2>d0)-.5)*1);
            curr_theta = curr_deg./180.*pi;
            theta= -1*curr_theta+self.theta_offset;
            
            normal = Q1 - Q2;
            v=null(normal);
            radius = curr_r;
            cart_coords=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
            cart_coords = cart_coords';
            
            %error('@@')
        end
        
        function cart_coords = cyl2cart_bat(self,cyl_coords); % batch conversion, added 210703
            d0 = self.d_0;
            Q1 = self.Q_1;
            Q2 = self.Q_2;
            
            len = cyl_obj.check_nby3(cyl_coords);
            curr_r = cyl_coords(:,1);
            curr_deg = self.angle_sign*cyl_coords(:,2);
            curr_z = cyl_coords(:,3);
 
            
            center = [self.k1*curr_z+self.b1,self.k2*curr_z+self.b2,self.k3*curr_z+self.b3];
            curr_theta = curr_deg./180.*pi;
            theta= -1.*curr_theta+self.theta_offset;
            normal = Q1 - Q2;
            v=null(normal);
            
            %v = repmat(v,[len,1]);
            
            %cart_coords=repmat(center',1,size(theta,2))+curr_r*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
            
            theta_ext = repmat(theta,[1,3]);
            cart_coords = center + repmat(curr_r,[1,3]).*(repmat(v(:,1)' ,[len,1]).*cos(theta_ext)+ repmat(v(:,2)' ,[len,1]).*sin(theta_ext) );
            %cart_coords=center+curr_r*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
        end
        
        function set.angle_sign(self,val) % added 210528
            if ~isnumeric(val)
                error('invalid input!')
            end
            if ~isscalar(val)
                error('invalid input!')
            end
            if ~ismember(val,[-1,1])
                error('angle_sign must be 1 or -1!')
            end
            self.angle_sign = val;
            
        end
        
        function out = get_x(self,cyl_coords) % these are for function handles
            cart_coords = cyl2cart(self,cyl_coords);
            out = cart_coords(1); 
        end
        function out = get_y(self,cyl_coords) % these are for function handles
            cart_coords = cyl2cart(self,cyl_coords);
            out = cart_coords(2); 
        end        
        function out = get_z(self,cyl_coords) % these are for function handles
            cart_coords = cyl2cart(self,cyl_coords);
            out = cart_coords(3); 
        end
        function show_cyl(self) 

            
            theta = linspace(-180,180,1000);
            theta = theta(:);
            
            
            cyl_pos = zeros(length(theta),3);
            cyl_pos(:,2) = theta;
            cyl_pos(:,1) = .5.*self.d_0;
            cyl_pos(:,3) = .5.*self.d_0;
            
            cart_pos = self.cyl2cart_bat(cyl_pos);
            

            
            %tmp_ref_pt = self.cart2cyl(self.ref_pt);
            tmp_ref_pt = self.cyl2cart([.5.*self.d_0,0,.5.*self.d_0]);
            
            
            new_pos = cart_pos;
            new_pos = [new_pos;self.Q_1;self.Q_2;tmp_ref_pt];
            
                        figure(1007)
                        view(3)
                        xlim([min(new_pos(:,1)),max(new_pos(:,1))]);
            ylim([min(new_pos(:,2)),max(new_pos(:,2))]);
            zlim([min(new_pos(:,3)),max(new_pos(:,3))]);
            hold on
            
            arrow3(self.Q_1,self.Q_2)
                        traj3D(cart_pos(:,1),cart_pos(:,2),cart_pos(:,3));
            sct1 = scatter3(tmp_ref_pt(1),tmp_ref_pt(2),tmp_ref_pt(3));
            
            sct1.Marker = 'o';
            sct1.LineWidth = 2;
            sct1.MarkerEdgeColor = [100 149 237]./256;
            sct1.MarkerFaceColor = [100 149 237]./256;
            sct1.SizeData = 10;

            %axis tight
            
            
        end
        
        function tf = is_rh(self) % is current cyl obj a right handed helix?
            new_p1 = self.cyl2cart([0,0,0]);
            new_p2 = self.cyl2cart([0,0,1]);
            new_q  = self.cyl2cart([1,0,0]);
            new_ref = self.cyl2cart([1,-90,0]);
            
            new_cyl = cyl_obj([0,0,0],new_p2 - new_p1,new_q-new_p1);
            
            temp_cyl = new_cyl.cart2cyl(new_ref-new_p1);
            
            if abs(temp_cyl(2)-90)<1
                tf = false;
            elseif abs(temp_cyl(2)+90)<1
                tf = true;
            else
                tf = true;
                warning('unexpected result')
            end
            
        end
        
    end
    
    methods(Static)
        function data_len = check_nby3(in) % check if input is an n by 3 array and return length
            if ~isnumeric(in)
                error('invalid input!')
            end
            sz = size(in);
            if length(sz)~=2
                error('invalid input!')
            end
            data_len = length(in(:,1));
            
        end
        function out = pt2line(P,Q,R) 
            
            
            len1 = cyl_obj.check_nby3(P);
            len2 = cyl_obj.check_nby3(Q);
            len3 = cyl_obj.check_nby3(R);
            if len1~=len2 || len2~=len3
                error('inputs inconsistent in length!')
            end

            %t =  sum((P - Q).*(R - P))./sum((P-Q).^2);
            t =  sum((P - Q).*(R - P),2)./repmat(sum((P-Q).^2,2), [1,3]);
            out = t.*(P-Q) + P;
            
            %error('!!')
        end
        
    end
    
    
    
end

% test commands

%% this returns a cyl object
% cyl0 = cyl_obj([32.4903   36.6040   19.6485],[33.2868   36.6679   20.6547],[33.1811   36.4856   20.5232]);
%% this returns cyl coordinates of [20,10,-50] in cartesian 
% cart2cyl(cyl0,[20,10,-50])
% ans = 40.8502  -55.3485  -63.6078
%% converting the cyl coords back to cart
% cyl2cart(cyl0,[ 40.8502  -55.3485  -63.6078])
% ans = 19.9996    9.9993  -49.9997

