classdef lf_obj < handle % linear fit object
   properties (SetAccess = protected)
       x_data
       y_data
   end
   properties (SetAccess = protected) % these should not be dependent, takes too much resources
       k
       b
       x_1
       y_1
   end
   
   properties % with set options
      % fit mode must be 1, 2 or 3
      % 1: best r_sqr
      % 2: best k (abs)
      % 3: y closest to zero
       fit_mode 
       
       N % how many x data points used for fitting?
   end
   methods
       function self = lf_obj(varargin)
               self.x_data = varargin{1};
               self.y_data = varargin{2};
           if ~isvector(self.x_data)||~isvector(self.y_data) ... 
                   ||length(self.x_data)~=length(self.y_data)||length(self.x_data)<3
              error(' invalid inputs!')
           else
               self.x_data = self.x_data(:);
               self.y_data = self.y_data(:);
           end           
           
           if nargin == 2
               self.N = length(self.x_data);
           elseif nargin == 3
               self.N = varargin{3};            
           else 
               error('must be 2 or 3 inputs')
           end

           
       end
       function update_params(self)
           if ~isempty(self.fit_mode) && ~isempty(self.N)
                  [self.x_data,idx] = sort(self.x_data);
                  self.y_data = self.y_data(idx);
               if self.fit_mode == 3
                  dist_to_0 = abs(self.y_data);
                  [~,dist_idx] = sort(dist_to_0); % modified 200815: these should be continuous
                  
                  if any( unique(diff(sort(dist_idx)))~=1)
                     warning(['when trying to locate ' num2str(self.N) ...
                         ' y values closest to 0 it was found that the indices were not continuous'])
                     if length(intersect(find(dist_idx> dist_idx(1)) , find(dist_idx< dist_idx(1)+self.N-1))) >= length(intersect(find(dist_idx> dist_idx(end)-self.N+1) , find(dist_idx< dist_idx(end))))
                        dist_idx = dist_idx(1):dist_idx(1)+self.N-1;
                     else
                         dist_idx = dist_idx(end)-self.N+1:dist_idx(end);
                     end
                      
                  end
                  
                  self.x_1 = self.x_data(dist_idx(1:self.N));
                  self.y_1 = self.y_data(dist_idx(1:self.N));
                  
                  [self.k,self.b,~] = lf_obj.linear_fit_single(self.x_1,self.y_1);
                  if ~(any(self.y_1>0)&&any(self.y_1<0))
                     warning('something is wrong'); 
                  end
                  
               else
                   temp_k = []; temp_b = [];temp_r2 = [];
                   for i = 1:length(self.x_data)+1-self.N
                       [temp_k(i),temp_b(i),temp_r2(i)] = lf_obj.linear_fit_single(self.x_data(i:i+self.N-1), ... 
                           self.y_data(i:i+self.N-1));
                       
                   end
                   if self.fit_mode == 1
                       temp_idx = find(temp_r2 == max(temp_r2));
                       temp_idx = temp_idx(1);
                       self.x_1 = self.x_data(temp_idx:temp_idx+self.N-1);
                       self.y_1 = self.y_data(temp_idx:temp_idx+self.N-1);
                       self.k = temp_k(temp_idx);
                       self.b = temp_b(temp_idx);
                   elseif self.fit_mode == 2
                       temp_idx = find(abs(temp_k) == max(abs(temp_k)));
                       temp_idx = temp_idx(1);
                       self.x_1 = self.x_data(temp_idx:temp_idx+self.N-1);
                       self.y_1 = self.y_data(temp_idx:temp_idx+self.N-1);
                       self.k = temp_k(temp_idx);
                       self.b = temp_b(temp_idx);                       
                   else
                      error('something is wrong!') 
                   end
               end
           end
       end
       function show_fitting(self)
           figure(1004)
           hold off
           sct1 = scatter(self.x_data,self.y_data);
           hold on
           sct1.Marker = 's';
            sct1.LineWidth = 1.2;
            sct1.SizeData = 12;
            sct1.MarkerEdgeColor = [100 149 237]./256;
            sct1.MarkerFaceColor = [100 149 237]./256;
          pp = plot(self.x_1,self.k.*self.x_1 + self.b);
          pp.Color = [255,140,0]./256;
          pp.LineWidth = 1.2;
          title([' k = ' num2str(self.k)])
       end
       function set.N(self,val)
          if ~isnumeric(val)|| ~isscalar(val)
              error('N must be a scalar!')
          elseif round(val)~=val
              error('N must be an integer!')
          elseif val>length(self.x_data)
              error('N must be <= length of x/y')
          else
              self.N = val;
              self.update_params
          end
       end
       function set.fit_mode(self,val)
           if ~ismember(val,[1,2,3])
              error('mode must be 1,2 or 3!')
           else
              self.fit_mode = val; 
              self.update_params;
           end
       end
       
   end
   methods(Static)
       function [k,b,r2] = linear_fit_single(x,y) 
            [P,S] = polyfit(x,y,1);

            k = P(1);
            b = P(2);
            r2 = 1 - (S.normr/norm(y - mean(y)))^2;
       end
   end
    
end