function F = cylFun(x,xp,yp,zp)

xi=[xp yp zp];
o=repmat([x(1) x(2) x(3)],length(xp),1);
a=[x(4) x(5) x(6)];
r=repmat(x(7),length(xp),1);
ri=xi-o-(xi*a')*a;

di=sqrt(ri(:,1).^2+ri(:,2).^2+ri(:,3).^2)-r;
F = sum(di.^2);