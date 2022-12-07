% get range and axis limits for input data
% x/y: x_1d, y_1d

% modified 210518: enabled struct input

% old calling:data = get_img_dims(x,y,x_step,y_step)

function data = get_img_dims(varargin)

if nargin == 1
    in = varargin{1};
    x = in.x_ctr;
    y = in.y_ctr;
    x_step = in.x_step;
    y_step = in.y_step;
    
elseif nargin == 4
   x = varargin{1};
   y = varargin{2};
   x_step = varargin{3};
   y_step = varargin{4};
    
else
    error('unknown error!')
end

data.x_lim = [min(x(:)) - .5*x_step,max(x(:)) + .5*x_step];
data.x_rg = data.x_lim(2)-data.x_lim(1);

data.y_lim = [min(y(:)) - .5*y_step,max(y(:)) + .5*y_step];
data.y_rg = data.y_lim(2)-data.y_lim(1);


end