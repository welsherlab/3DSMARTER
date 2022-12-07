function draw_dashed_circle(x,y,r ,N,lw,clr)
th = linspace(0,2*pi,N*2+1);
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;

for i = 1:N
   pp = plot(xunit(i*2-1:i*2) ,yunit(i*2-1:i*2));
    pp.LineWidth = lw;
    pp.Color = clr;
    
end

end