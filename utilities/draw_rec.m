% draw a rectangle frame, returns plot object

function p = draw_rec(x_rg,y_rg,lw,clr)
hold on

x_rg = x_rg(:);
y_rg = y_rg(:);

x1 = min(x_rg);
x2 = max(x_rg);

y1 = min(y_rg);
y2 = max(y_rg);

p = plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1]);
p.LineWidth = lw;
p.Color = clr;

end