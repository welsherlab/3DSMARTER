function [c,ceq] = cylCon(x)
c=[];
ceq=[];
ceq(1) = dot([x(1) x(2) x(3)],[x(4) x(5) x(6)]);
ceq(2)= norm([x(4) x(5) x(6)])-1;
