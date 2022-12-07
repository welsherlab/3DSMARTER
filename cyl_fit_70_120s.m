clear;cab

load(fullfile(pwd, '\data', '200304 TR003 filopodia 1usec 70_120s.mat'));
[params ,fval ,hessian ,xx, yy, zz]=parsedCylFit(xOut(1:1000:end),yOut(1:1000:end),zOut(1:1000:end));

% obtain center axis and build cylinder object

Q1 = [mean(xx(1,:)),mean(yy(1,:)),mean(zz(1,:))];
Q2 = [mean(xx(2,:)),mean(yy(2,:)),mean(zz(2,:))];

P = [10,20,30]; % arbitraty ref point (must be out of line Q1-Q2)

cyl = rcyl(Q1, Q2, P);

% downsample coordinate data, convert to coordinates in cylindrical system

x = xOut(1:1e3:end);
y = yOut(1:1e3:end);
z = zOut(1:1e3:end);
%%
cyl_pos = cyl.cart2cyl([x(:), y(:), z(:)]);

% plot r, theta, h positions
t = [1:length(cyl_pos)]./1e3;

figure;
plot(t, cyl_pos(:,1));
xlabel('t [sec]')
ylabel('r [\mum]')

figure;
plot(t, unwrap(cyl_pos(:,2)));
xlabel('t [sec]')
ylabel('theta [бу]')

figure;
plot(t, cyl_pos(:,3));
xlabel('t [sec]')
ylabel('h [\mum]')
