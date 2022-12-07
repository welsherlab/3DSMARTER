function [Dmap,VmapX,VmapY]=cylSurfProps
%Analyze the diffusion and force properties of unwrapped filopodial
%coordinates
%%210529 KDW

%% Load data
[fname,pname]=uigetfile('*.mat','Open "unwrapped 2d traj 590-650s.mat"');
load([pname fname])
%% Grid data
%Decimate to 1 msec
xCyl=decimate(x,1e3);
yCyl=decimate(y,1e3);

%Discretize xy values
[xGridPos,xb]=discretize(xCyl,60);
[yGridPos,yb]=discretize(yCyl,60);
bincentersX=mean([xb(1:end-1);xb(2:end)]);
bincentersY=mean([yb(1:end-1);yb(2:end)]);
%% Get diffusion map
noise=0.006;
[density,Dmap,D_global]=diffusionMap(x,y,noise,xGridPos,yGridPos,bincentersX,bincentersY);

%% Get force map along y
[densityV,VmapY]=forceMap(x,y,D_global,noise,xGridPos,yGridPos,bincentersX,bincentersY);

%% Get force map along x
[densityX,VmapX]=forceMapX(x,y,r,D_global,noise,xGridPos,yGridPos,bincentersX,bincentersY);

%% Get spring constant and potential
cutoff=25;
getScalarPotential(VmapY,VmapX,density,D_global,cutoff,bincentersX,bincentersY)

function [density,Dmap,D_global]=diffusionMap(x,y,noise,xGridPos,yGridPos,bincentersX,bincentersY)
%Generate force map for cylinder surface diffusion
%210529 KDW

%set up range of Ds to test
D_basis=0.01:0.01:.5;
D=repmat(D_basis,max(xGridPos),1,max(yGridPos));

%% Initialize figures
Dfig=figure('Name','D map');
%%
%recursive Bayesian inference to find the diffusion coefficient
%Calculate steps in x and y, try to do this for multiple values of tau
for tau=1e-3:1e-3:1e-2
    %Get steps along y only. This is to avoid unrapping issues along the x
    %direction
    ystep=diff(decimate(y,round(tau/1e-6)));
    
    %calculated sigma_squared
    sigma_squared=2*D*tau+noise^2;
    
    %Initialize prior only for first time step
    if round(tau,4)==1e-3
        prior=ones(size(D));
        prior=prior./sum(prior,2);
        lPr=log(prior);
        sum(exp(lPr),2);
    end
    
    %Loop through displacements
    for i=1:length(ystep)
        
        %Find the current grid position
        xi=xGridPos(i);
        yi=yGridPos(i);
        
        %Generate a mask to place the likelihood of the current step into
        %the gridded map
        mask=zeros(size(lPr));
        mask(xi,:,yi)=ones(1,size(lPr,2),1);
        
        %calculate current likelihood (log)
        ll=-0.5*log(sigma_squared)-ystep(i).^2./2./sigma_squared;
        
        %Calculate posterior
        lPo=mask.*ll+lPr;
        
        %normalize posterior
        lPo=lPo-max(lPo,2);
        lPo=lPo-log(sum(exp(lPo),2));
        
        %Posterior becomes prior
        lPr=lPo;
        
        %Normalization check
        checknorm=sum(exp(lPo),2);
        if ~any(round(checknorm,3)==1)
            disp('norm error')
            break
        end
    end
    figure(Dfig)
    [~,in]=max(lPo,[],2);
    imagesc(bincentersY,bincentersX,squeeze(D_basis(in)))
    pause(0.01)
    drawnow
end
%%
%Generate D map
Dmap=squeeze(D_basis(in));

%Remove undersampled regions
density=accumarray([xGridPos yGridPos],1);
Dmap(density<10)=NaN;
figure(Dfig)
imagesc(bincentersY,bincentersX,Dmap)
D_global=nanmean(Dmap,'all');

figure('Name','Density')
imagesc(bincentersY,bincentersX,density)

function [density,Vmap]=forceMap(x,y,D_global,noise,xGridPos,yGridPos,bincentersX,bincentersY)
%Generate force map for cylinder surface diffusion
%210529 KDW
%set up range of Vs to test
V_basis=(-10:.1:10);
V=repmat(V_basis,max(xGridPos),1,max(yGridPos));

%%
%recursive Bayesian inference to find the drift velocity
%tau = 1 ms
tau=1e-3;
%Initialize prior
prior=ones(size(V));
prior=prior./sum(prior,2);
lPr=log(prior);
sum(exp(lPr),2);

%Get steps along y only. This is to avoid unrapping issues along the x
%direction
ystep=diff(decimate(y,round(tau/1e-6)));

%calculated sigma_squared
sigma_squared=2*D_global*tau+noise^2;

%Loop through displacements
for i=1:length(ystep)
    
    %Find the current grid position
    xi=xGridPos(i);
    yi=yGridPos(i);
    
    %Generate a mask to place the likelihood of the current step into
    %the gridded map
    mask=zeros(size(lPr));
    mask(xi,:,yi)=ones(1,size(lPr,2),1);
    
    %calculate current likelihood (log)
    ll=-0.5*log(sigma_squared)-(ystep(i)-V.*tau).^2./2./sigma_squared;
    
    %Calculate posterior
    lPo=mask.*ll+lPr;
%     plot(V_basis,lPo(xi,:,yi))
%     pause(0.01)
%     drawnow
    %normalize posterior
    lPo=lPo-max(lPo,2);
    lPo=lPo-log(sum(exp(lPo),2));
    
    %Posterior becomes prior
    lPr=lPo;
    
    %Normalization check
    checknorm=sum(exp(lPo),2);
    if ~any(round(checknorm,3)==1)
        disp('norm error')
        break
    end
end
%%
%Generate V map
[~,in]=max(lPo,[],2);
Vmap=squeeze(V_basis(in));

%Calculate density
density=accumarray([xGridPos yGridPos],1);

function [density,Vmap]=forceMapX(x,y,r,D_global,noise,xGridPos,yGridPos,bincentersX,bincentersY)
%Generate force map for cylinder surface diffusion
%210529 KDW

%set up range of Vs to test
V_basis=(-10:.1:10);
V=repmat(V_basis,max(xGridPos),1,max(yGridPos));
%%
%recursive Bayesian inference to find the drift velocity
%tau = 1 ms
tau=1e-3;
%Initialize prior
prior=ones(size(V));
prior=prior./sum(prior,2);
lPr=log(prior);
sum(exp(lPr),2);

%Get steps along x only.
%Get angle
theta=unwrap(x./mean(r));
xstep=diff(decimate(mean(r)*theta,round(tau/1e-6)));
%There are a significant number of
%data points where theta = +/- 2pi. For force mapping, we should simply
%skip these points
% xcut=x;
% ycut=y;
% xcut([abs(diff(theta))>(1.5*pi);false])=NaN;
% ycut([abs(diff(theta))>(1.5*pi);false])=NaN;
% Cx=breakNaNData(xcut);
% Cy=breakNaNData(ycut);
% density=zeros(60,60);

%% Loop through displacements
%calculated sigma_squared
sigma_squared=2*D_global*tau+noise^2;
% for j=1:length(Cx)
%     %only use data with more than 5 ms
%     if length(Cx{j})>5e3
%         tempX=Cx{j};
%         tempY=Cy{j};
%         xstep=diff(decimate(tempX,round(tau/1e-6)));
for i=1:length(xstep)
    %Find the current grid position
    xi=xGridPos(i);
    yi=yGridPos(i);
    %             %Find the current grid position
    %             [~,xi]=min(abs(bincentersX-tempX(i)));
    %             [~,yi]=min(abs(bincentersY-tempY(i)));
    %             density(xi,yi)=density(xi,yi)+1;
    %Generate a mask to place the likelihood of the current step into
    %the gridded map
    mask=zeros(size(lPr));
    mask(xi,:,yi)=ones(1,size(lPr,2),1);
    
    %calculate current likelihood (log)
    ll=-0.5*log(sigma_squared)-(xstep(i)-V.*tau).^2./2./sigma_squared;
    
    %Calculate posterior
    lPo=mask.*ll+lPr;
    
    %normalize posterior
    lPo=lPo-max(lPo,2);
    lPo=lPo-log(sum(exp(lPo),2));
    %                 plot(V_basis,lPo(xi,:,yi))
    %     pause(0.01)
    %     drawnow
    %Posterior becomes prior
    lPr=lPo;
    
    %Normalization check
    checknorm=sum(exp(lPo),2);
    if ~any(round(checknorm,3)==1)
        disp('norm error')
        break
    end
end
%         imagesc(density)
%         axis image
%         drawnow
%     end
% end
%%
%Generate V map
[~,in]=max(lPo,[],2);
Vmap=squeeze(V_basis(in));

%Calculate density
density=accumarray([xGridPos yGridPos],1);

function getScalarPotential(Vmap,VmapX,density,D_global,cutoff,bincentersX,bincentersY)
%Find drag coefficient
kb=1.381e-23;
T=310.15;
gamma=kb*T/D_global/1e-12;
%Generate positive quivers
qPos=Vmap*gamma*1e-6*1e12;
qPos(qPos<0)=NaN;

qPos(density<cutoff)=NaN;
%Generate negative quivers
qNeg=Vmap*gamma*1e-6*1e12;
qNeg(qPos>0)=NaN;
qNeg(density<cutoff)=NaN;
[qX,qY]=meshgrid(bincentersY,bincentersX);
figure('Name','Force Map Y')
imagesc(bincentersY,bincentersX,density)
hold on
quiver([qX(:);min(bincentersY)],[qY(:);max(bincentersX)],[qPos(:);.1],[zeros(size(Vmap(:)));0],'Color','g','LineWidth',1)
quiver([qX(:);max(bincentersY)],[qY(:);max(bincentersX)],[qNeg(:);-.1],[zeros(size(Vmap(:)));0],'Color','r','LineWidth',1)
hold off
axis ij
axis image
% Now X

%Generate positive quivers
qPosX=VmapX*gamma*1e-6*1e12;
qPosX(qPosX<0)=NaN;
qPosX(density<cutoff)=NaN;
%Generate negative quivers
qNegX=VmapX*gamma*1e-6*1e12;
qNegX(qPosX>0)=NaN;
qNegX(density<cutoff)=NaN;

figure('Name','Force Map X')
imagesc(bincentersY,bincentersX,density)
hold on
quiver([qX(:);min(bincentersY)],[qY(:);max(bincentersX)],[zeros(size(Vmap(:)));0],[qPosX(:);.1],'Color','g','LineWidth',1)
quiver([qX(:);max(bincentersY)],[qY(:);max(bincentersX)],[zeros(size(Vmap(:)));0],[qNegX(:);-.1],'Color','r','LineWidth',1)
hold off
axis ij
axis image

% Now both
%Generate horizontal quivers
qHorz=Vmap*gamma*1e-6*1e12;
qHorz(density<cutoff)=NaN;
%Generate vertical quivers
qVert=VmapX*gamma*1e-6*1e12;
qVert(density<cutoff)=NaN;
% figure('Name','V map XY White')
% imagesc(bincentersY,bincentersX,density)
% hold on
% quiver([qX(:);min(bincentersY)],[qY(:);max(bincentersX)],[qHorz(:);.1],[qVert(:);.1],'Color','w','LineWidth',1)
% hold off
% axis ij
% axis image
% colormap plasma

% Determine angles for colormapping
%Generate horizontal quivers
qHorz=Vmap(:)*gamma*1e-6*1e12;
qHorz(density<cutoff)=NaN;
%Generate vertical quivers
qVert=VmapX(:)*gamma*1e-6*1e12;
qVert(density<cutoff)=NaN;
phi=atan2(qHorz,qVert);

n=96;   %cmap elements
cmapstep=range(phi)/n;
mincmap=min(phi);
maxcmap=max(phi);
mincmap=-pi;
maxcmap=pi;

%build a circular colormap
startColor=[1 0 0]
secondColor=[0 1 0];
thirdColor=[0 0 1];
finalColor=startColor;

mFirstRamp=(secondColor-startColor)/round(n/3);
firstRamp=repmat((1:round(n/3))',1,3).*repmat(mFirstRamp,round(n/3),1)+startColor;

mSecondRamp=(thirdColor-secondColor)/round(n/3);
secondRamp=repmat((1:round(n/3))',1,3).*repmat(mSecondRamp,round(n/3),1)+secondColor;

mThirdRamp=(finalColor-thirdColor)/round(n/3);
thirdRamp=repmat((1:round(n/3))',1,3).*repmat(mThirdRamp,round(n/3),1)+thirdColor;

basemap=[firstRamp;secondRamp;thirdRamp];
cmap=vertcat(basemap);
vectorScale=.1;

hVmap=figure('Name','Force Map XY');
imagesc(bincentersY,bincentersX,density)
hold on
for i=1:n
    %get current color
    currentMin=mincmap+(i-1)*cmapstep;
    currentMax=mincmap+(i)*cmapstep;
    mask=(phi>=currentMin)&(phi<currentMax);
    quiver(qX(mask),qY(mask),qHorz(mask)*vectorScale,qVert(mask)*vectorScale,0,'Color',cmap(i,:),'LineWidth',0.5)
    if ~mod(i-1,6)
        scaleArrowLength=0.5*vectorScale;  %50 fN scale arrow
        arrowHeadScale=.01;
    quiver([min(bincentersY)+.1 0],[max(bincentersX)-.225 0],[scaleArrowLength*sin(mean([currentMin currentMax])) arrowHeadScale],[scaleArrowLength*cos(mean([currentMin currentMax])) arrowHeadScale],0,'Color',cmap(i,:),'LineWidth',0.5)
    end
end

hold off
colormap gray
axis image
%100 nm scale bar
h=line([min(bincentersY)+.05 min(bincentersY)+.15],[max(bincentersX)-.16 max(bincentersX)-.16],'Color','w','LineWidth',2);
axis([0 1 -.5 .5])
h=rectangle('Position',[0.035 0.1 .35 .35],'EdgeColor','w','LineWidth',2,'LineStyle','--');
axis off
fh = copyobj(hVmap, 0);
h=findobj(gca,'type','rectangle');
delete(h)
axis([0.035 .385 .1 .45])
% Line integrals!
%Generate horizontal quivers
qHorz=Vmap*gamma*1e-6*1e12;

%Generate vertical quivers
qVert=VmapX*gamma*1e-6*1e12;

dx=mean(diff(bincentersY));
dy=mean(diff(bincentersX));
%figure('Name','Horizontal Potential');
%imagesc(bincentersY,bincentersX,qHorz)
%qHorz - x, but goes with bincentersY, rows
%qVert - y, but goes with bincentersX, columns

qHorz(density<cutoff)=0;
%Integrate qHorz first
vHorz=-cumtrapz(dx,qHorz,2);
vHorz(density<cutoff)=0;
%imagesc(bincentersY,bincentersX,vHorz)

%figure('Name','Vertical Potential')
qVert(density<cutoff)=0;
vVert=-cumtrapz(dy,qVert,1);
vVert(density<cutoff)=0;
%imagesc(bincentersY,bincentersX,vVert)

%need a reference point
% figure('Name','Mean Potential')
normMesh=(vHorz+vVert)/2;
% surf(bincentersY,bincentersX,normMesh*1e-18/kb/T)
% axis square

[xMesh,yMesh]=ndgrid(bincentersX,bincentersY);
%% Fit to surface mesh
[xData, yData, zData] = prepareSurfaceData( xMesh, yMesh, normMesh );

% Set up fittype and options.
ft = 'thinplateinterp';

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );

% Plot fit with data.
figure( 'Name', 'Potential Surface Fit' );
h = plot( fitresult, [xData, yData], zData );
legend( h, 'Potential Surface Fit', 'Mean Potential', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'xMesh', 'Interpreter', 'none' );
ylabel( 'yMesh', 'Interpreter', 'none' );
zlabel( 'normMesh', 'Interpreter', 'none' );
grid on
view( 55.9, 42.6 );
lineHandles=findobj(gca,'type','line');
lineHandles.MarkerSize=3;
