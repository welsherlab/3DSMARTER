function [params ,fval ,hessian ,xOut, yOut, zOut]=parsedCylFit(xp,yp,zp)


shortXP=xp(:);
shortYP=yp(:);
shortZP=zp(:);
xo=[mean(shortXP) mean(shortYP) mean(shortZP) mean(shortXP) mean(shortYP) mean(shortZP) .3];
LB=[-100 -100 -100 -100 -100 -100 0];
UB=[100 100 100 100 100 100 10];
options = optimset('Display','iter','Algorithm','active-set','PlotFcns',@optimplotx);
exitflag=0;
while exitflag==0
    [xo,fval,exitflag]=fmincon(@(x) cylFun(x,shortXP,shortYP,shortZP),...
        xo,[],[],[],[],LB,UB,@(x) cylConFixR(x,xo(7)),options);
end
exitflag=0;
while exitflag==0
    [xo,fval,exitflag,~,~,~,hessian]=fmincon(@(x) cylFun(x,shortXP,shortYP,shortZP),...
        xo,[],[],[],[],LB,UB,@cylCon,options);
end
x=xo;
params=x';
o=x(1:3);
a=x(4:6);
r=x(7);
[xOut,yOut,zOut]=drawCyl(o,a,r,shortXP,shortYP,shortZP,1);
