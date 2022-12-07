function data = ct_3(varargin)
% count spots in 2d or 3d cartesian system in specified grid

if nargin == 4 % x,y,x_grid(x_step),y_grid(y_step)
    x = check_input(varargin{1});
    y = check_input(varargin{2});
    
    [x_grid,x_step] = check_input2(varargin{3},x);
    
    [y_grid,y_step] = check_input2(varargin{4},y);
    
    z = zeros(length(x),1);
    z_grid = [-.5,.5];
    z_step = 1;
    z_l = 1;
    
elseif nargin == 6 % x,y,z,x_grid(x_step),y_grid(y_step),z_grid(z_step)
    x = check_input(varargin{1});
    y = check_input(varargin{2});
    z = check_input(varargin{3});
    
    [x_grid,x_step] = check_input2(varargin{4},x);
    [y_grid,y_step] = check_input2(varargin{5},y);
    [z_grid,z_step] = check_input2(varargin{6},z);
    
else
    
    error('insert 4 or 6 inputs!')
    
end

data.x_step = x_step;
data.y_step = y_step;
data.z_step = z_step;

x_l = length(x_grid)-1;
y_l = length(y_grid)-1;
z_l = length(z_grid)-1;

% evaluate x,y,z length

if length(x)~=length(y) 
error('data length does not match!')
end

if length(x)~=length(z) 
error('data length does not match!')
end

% first, remove out of range data points, show warning if there are such
% points

bad_x = union(find(x<=x_grid(1)),find(x>=x_grid(end)));
bad_y = union(find(y<=y_grid(1)),find(y>=y_grid(end)));
bad_z = union(find(z<=z_grid(1)),find(z>=z_grid(end)));

bad_idx = union(union(bad_x,bad_y),bad_z);
if ~isempty(bad_idx)
    warning([num2str(length(bad_idx)) ' out of ' num2str(length(x)) ' points are out of range!'])
    x(bad_idx) = [];
    y(bad_idx) = [];
    z(bad_idx) = [];
end

x_idx = ceil((x-x_grid(1))/x_step);
y_idx = ceil((y-y_grid(1))/y_step);
z_idx = ceil((z-z_grid(1))/z_step);  

idx3 = [x_idx(:),y_idx(:),z_idx(:)];

%save('idx3 probe.mat','idx3')

data.idx3 = idx3;

x_ctr = movmean2(x_grid);
y_ctr = movmean2(y_grid);
z_ctr = movmean2(z_grid);

data.xg = x_ctr;
data.yg = y_ctr;
data.zg = z_ctr;

%[ineq_row,row_ct] = ct_row(idx3);
%[length(x_ctr),length(y_ctr),length(z_ctr)]
[ineq_row,row_ct] = ct_idx3(idx3,[length(x_ctr),length(y_ctr),length(z_ctr)]);

data.rows = ineq_row;
data.row_ct = row_ct;

temp = zeros(length(row_ct),1);

data.x_ctr = temp;
data.y_ctr = temp;
data.z_ctr = temp;

data.int = temp;

for i = 1:length(temp)
    data.x_ctr(i) = x_ctr(ineq_row(i,1));
    data.y_ctr(i) = y_ctr(ineq_row(i,2));
    data.z_ctr(i) = z_ctr(ineq_row(i,3));    
    data.int(i)    = row_ct(i);
end

% I don't understand why these two vars need to be switched, but it seems
% this is the right way
[data.x_ctr_all,data.y_ctr_all,data.z_ctr_all] = meshgrid(y_ctr,x_ctr,z_ctr);

data.int_all = zeros(size(data.x_ctr_all));

for i = 1:length(temp)
    data.int_all(ineq_row(i,1),ineq_row(i,2),ineq_row(i,3)) = row_ct(i);
end



end

function out = check_input(in) % input must be a vector
    if isscalar(in)
        error('input must be a vector!')
    end
    if ~isvector(in)
        error('invalid input!')
    else
        out = in(:);
    end
end

function [out,out_step ]= check_input2(in1,in2) % in1: input as x_grid or x_step in2: input x points 
    if isscalar(in1)
        out = ezgrid(in2,in1);
        out_step = in1;
    else
        out = check_input(in1);
        out_step = mean(diff(in1));
    end
end