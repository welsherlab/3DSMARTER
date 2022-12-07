function [f,v] = get_patch_data_2d(x_grid_1d,y_grid_1d,x_step,y_step)

L = length(x_grid_1d);

v = zeros(4*length(x_grid_1d),2)-10;

v([1:L].*4-3,1) = x_grid_1d - .5.*x_step; 
v([1:L].*4-3,2) = y_grid_1d + .5.*y_step;

v([1:L].*4-2,1) = x_grid_1d + .5.*x_step; 
v([1:L].*4-2,2) = y_grid_1d + .5.*y_step;

v([1:L].*4-1,1) = x_grid_1d + .5.*x_step; 
v([1:L].*4-1,2) = y_grid_1d - .5.*y_step;

v([1:L].*4-0,1) = x_grid_1d - .5.*x_step; 
v([1:L].*4-0,2) = y_grid_1d - .5.*y_step;

f = zeros(length(x_grid_1d),4);
for i = 1:4
f(:,i) = [([0:L - 1]).*4+i ]';
end