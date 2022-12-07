function traj3D_custom(xx,yy,zz,LW)
%figure;
len=length(xx);
cjet=colormap(jet(round(len/33)));

cjlen=length(cjet);
hold on;
for i=1:cjlen-1
    seg=(floor((i-1)*len/cjlen)+1):floor((i+1)*len/cjlen);
    p1=plot3(xx(seg),yy(seg),zz(seg));
    %%%%%
%     drawnow;
%     pause(0.1);
%     M(i) = getframe;
%     %%%%%%%%
    set(p1,'Color',cjet(i,:),'LineWidth',LW);
% set(p1,'Color',cjet(i,:),'LineWidth',2);
end
% movie 
%view(3);
%adjustPlotRange(gcf);