%v2: only dens plot, removed some data out of range, use circle rather than
% crossmark as bar
clear;cab

main_dir = 'Z:\Manuscripts\metal particle tracking\1 usec data\200304 TR003 filopodia\figures';

t0 = 400; % global t0
r_rg = [1,2]; % radius range to plot range of interest

% import MHz bayesian coordinates and cylindrical fitting data
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


point0 = [xOut(1),yOut(1),zOut(1)];

% get cylinder object
cyl = cyl_obj(Q1,Q2,point0);


cyl_pos = cyl.cart2cyl_bat(pos_1d);

cab
r_1 = mean(cyl_pos(:,1)); 

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

pos_1d = [x,y,z];
pos_1d(1,:) = [];
Q1 = [mean(xx(1,:)),mean(yy(1,:)),mean(zz(1,:))];
Q2 = [mean(xx(2,:)),mean(yy(2,:)),mean(zz(2,:))];


cyl = cyl_obj(Q1,Q2,point0);

cyl_pos = cyl.cart2cyl_bat(pos_1d);
cab
r_1 = mean(cyl_pos(:,1)); 
sim_pos = [];

% obtain 2d coordinates from top view

sim_pos(:,1) = cyl_pos(:,1).*cos(cyl_pos(:,2)./180.*pi);
sim_pos(:,2) = cyl_pos(:,1).*sin(cyl_pos(:,2)./180.*pi);

step_size = .01;


remove_idx = union(union(find(sim_pos(:,1)>.3),find(sim_pos(:,1)<-.3)),union(find(sim_pos(:,2)>.3),find(sim_pos(:,2)<-.3)));
sim_pos(remove_idx,:) = [];

xg =ezgrid(sim_pos(:,1),step_size);
yg = ezgrid(sim_pos(:,2),step_size);

% counting spots in designated grid
data = ct_3(sim_pos(:,1),sim_pos(:,2),xg ,yg);

[f,c,v] =  get_2d_pdata_from_ct3(data);



plt_h = 1.02;

clrs = Paired(12);
LW = .5;

plt = Plot([0,1],[0,1]);
plt.LineStyle = 'none';

hold on


axis equal
axis manual

img_dims = get_img_dims([-.3,.3],[-.3,.3],0,0);

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

fig = gcf;
fig.Units = 'pixels';
fig_pos = fig.Position; % [left bottom width height]
fig_pos(3) = fig_pos(3) + 100;
fig.Position = fig_pos;

p = draw_rec(plt.XLim,plt.YLim,LW,'k');

draw_dashed_circle(0,0,r_1 ,20,LW,clrs(8,:))

% add scale bar
pp = plot([-.05,.05],[0,0]);
pp.LineWidth = 1;
pp.Color = clrs(6,:);

pp = plot([0,0],[-.05,.05]);
pp.LineWidth = 1;
pp.Color = clrs(6,:);

