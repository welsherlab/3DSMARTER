% project coordiates to cylindrical surface and show density map

clear;cab

%%

t0 = 400; % global t0
r_rg = [1,2]; % radius range to plot range of interest

% load MHz bayesian position data and fitting parameters of cylindrical
% center axis
load('Z:\Manuscripts\metal particle tracking\1 usec data\200304 TR003 filopodia\200304 TR003 filopodia 1usec 400_460s corrected.mat')
load('Z:\Manuscripts\metal particle tracking\1 usec data\200304 TR003 filopodia\cyl data 400-460s.mat')

max_idx = length(xOut(1:1:end));
t_end = length(xOut)/1e6;

cab

clrs = Paired(12);

LW = 1;
f_sz = 8;

x = xOut(1:1:end);
y = yOut(1:1:end);
z = zOut(1:1:end);

plt = Plot([0,1],[0,1]);

hold on

plt.LineStyle = 'none';

pos_1d = [xOut(1:1:end),yOut(1:1:end),zOut(1:1:end)];
pos_1d(1,:) = [];
Q1 = [mean(xx(1,:)),mean(yy(1,:)),mean(zz(1,:))];
Q2 = [mean(xx(2,:)),mean(yy(2,:)),mean(zz(2,:))];

% pick an arbitrary reference point
point0 = [xOut(1),yOut(1),zOut(1)]; 

cyl = cyl_obj(Q1,Q2,point0);


cyl_pos = cyl.cart2cyl_bat(pos_1d);

cab
r_1 = mean(cyl_pos(:,1)); % average radius

%% trajectory 2d along surface

cab
sb_ctr = [ .35,.01];
sb_len = .1;

x = xOut(1:1:end);
y = yOut(1:1:end);
z = zOut(1:1:end);

plt = Plot([0,1],[0,1]);

hold on

plt.LineStyle = 'none';
cab

pos_1d = [x,y,z];
pos_1d(1,:) = [];
Q1 = [mean(xx(1,:)),mean(yy(1,:)),mean(zz(1,:))];
Q2 = [mean(xx(2,:)),mean(yy(2,:)),mean(zz(2,:))];

point0 = [xOut(1),yOut(1),zOut(1)];

cyl = cyl_obj(Q1,Q2,point0);

cyl_pos = cyl.cart2cyl_bat(pos_1d);

sim_pos(:,1) = 1.*pi.*r_1.*cyl_pos(:,2)./180;
sim_pos(:,2) = cyl_pos(:,3);

step_size = .01;

xg =ezgrid(sim_pos(:,1),step_size);
yg = ezgrid(sim_pos(:,2),step_size);

data = ct_3(sim_pos(:,1),sim_pos(:,2),xg ,yg);
[f,c,v] =  get_2d_pdata_from_ct3(data);

plt_h = 3.4;

clrs = Paired(12);
LW = 1;

plt = Plot([0,1],[0,1]);
plt.LineStyle = 'none';

hold on

img_dims = get_img_dims(data.x_ctr,data.y_ctr,step_size,step_size); %% this should be removed in sub plots

axis equal
axis manual

box_dim = [ img_dims.x_rg/img_dims.y_rg*plt_h,plt_h];

plt.BoxDim = box_dim;
plt.XLim = img_dims.x_lim;
plt.YLim = img_dims.y_lim;

axis off

ax = gca;
ax.Clipping = 'off';

c_map = plasma(max(c));

p_c = zeros(length(c),3);
for i = 1:length(c(:,1))
    p_c(i,:) = c_map(c(i),:);
    
end


curr_pt = patch('Faces',f,'Vertices',v,'FaceColor','flat','FaceVertexCData',p_c,'FaceAlpha',1,'LineStyle','none');
hold on



min_z = min(img_dims.y_lim);
max_z = max(img_dims.y_lim);



hold on

fig = gcf;
fig.Units = 'pixels';
fig_pos = fig.Position; % [left bottom width height]
fig_pos(3) = fig_pos(3) + 100;
fig.Position = fig_pos;

p = draw_rec(plt.XLim,plt.YLim,LW,'k');

