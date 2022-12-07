function [f,c,v] =  get_2d_pdata_from_ct3(data);

x_step = data.x_step;
y_step = data.y_step;

x = data.x_ctr;
y = data.y_ctr;
z = data.z_ctr;
c = data.int;

remove_idx = find(z>0);

x(remove_idx) = [];
y(remove_idx) = [];
c(remove_idx) = [];

[f,v] = get_patch_data_2d(x,y,x_step,y_step);

end