function [xOut,yOut,zOut]=drawCyl(o,a,r,xp,yp,zp,drawOpt)
projA=([xp yp zp]*a')*a;
[minVal,minIndex]=min([xp yp zp]*a');
[maxVal,maxIndex]=max([xp yp zp]*a');
[cx,cy,cz] = cylinder(1,2000);
cx=r.*cx;
cy=r.*cy;
cz=(maxVal-minVal).*cz;
size0=size(cx);
figure(1), clf;%surf(cx,cy,cy);axis equal;
%--- define Euler z-rotation matrix Rz
Rz=[o/norm(o);cross(a,o/norm(o));a];
P=[cx(:) cy(:) cz(:)]*Rz;  %rotate each point on surface
X=reshape(P(:,1),size0);%transform surface vertices back into matrix
Y=reshape(P(:,2),size0);
Z=reshape(P(:,3),size0);
clear P;

xOut=X+projA(minIndex,1)+o(1);
yOut=Y+projA(minIndex,2)+o(2);
zOut=Z+projA(minIndex,3)+o(3);
if drawOpt
    hold on;surf(xOut,yOut,zOut);axis equal;%light
    traj3D_custom(xp,yp,zp,1);
    shading interp; %lighting phong;axis tight
end
%note that rotation is performed through the origin