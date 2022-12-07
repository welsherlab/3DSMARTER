function [xOut,yOut,zOut]=bayesProcessing(ktPos,TAGPos,n,xStage,yStage,zStage)
%% Apply bayes filter
% Generate coords

%Diffusion params
D=.2;
tau=1e-6;

% Generate coords
%Generate XY coordinates
dx=.2e-3;
sigma=sqrt(2*D*tau);
noOfSigmas=150; %Number of std to include in localization grid
coords=gpuArray((-noOfSigmas*sigma):sigma:((noOfSigmas-1)*sigma));
[x,y]=meshgrid(coords);
z=gpuArray((-3*noOfSigmas*sigma):sigma:((noOfSigmas-1)*3*sigma));
diff=exp(-x.^2/4/D/tau).*exp(-y.^2/4/D/tau);    %Diffusion kernel
diffZ=exp(-z.^2/4/D/tau);    %Diffusion kernel Z
diff=diff((noOfSigmas-1):(noOfSigmas+3),(noOfSigmas-1):(noOfSigmas+3)); %cropped diffusion (two standard devations)
diffZ=diffZ((3*noOfSigmas-1):(3*noOfSigmas+3));
%%
%Initialize Prior
Pr=ones(size(x));
Pr=Pr/sum(Pr(:));
lPr=gpuArray(log(Pr));

%Initilize Z Prior
PrZ=ones(size(z));
PrZ=PrZ/sum(PrZ(:));
lPrZ=gpuArray(log(PrZ));

%Brightness and background
s=2e6;
b=2e6;

%Initial position
mu=0;

%Size of PSF
w=0.2;
wZ=0.5;

%Scan size
d=0.250;
xkt=gpuArray(d*[-1,1,1,-1]);
ykt=gpuArray(d*[1,1,-1,-1]);

%Loop through photons
index=1;
emission=0;

%Generate different possible "gammas" ahead of time. There are only 25, it
%does not need to be calculated each time through

lP_n_x=genLL(x,y,xkt,ykt,s,b,w,tau,n);
lP_n_z=genLLZ(z,TAGPos,s,b,wZ,tau,n);
%%
hWait=waitbar(0,'processing');
nData=1e5;
mleIndex=gpuArray(ones(nData,1));
mleIndexZ=gpuArray(ones(nData,1));
dataOut=[];
dataOutZ=[];
nit=100;
tic
for i=1:round(length(zStage)/nData)
    startIndex=(i-1)*nData+1;
    %define ktPos subset
    %define n subset
    [mleIndex,lPo]=bayesQuickCalc(lPr,lP_n_x,ktPos,n,diff,mleIndex,startIndex,nData);
    lPr=lPo;
    [mleIndexZ,lPoZ]=bayesQuickCalcZ(lPrZ,lP_n_z,TAGPos,n,diffZ,mleIndexZ,startIndex,nData);
    lPrZ=lPoZ;
    dataOut=vertcat(dataOut,gather(mleIndex));
    dataOutZ=vertcat(dataOutZ,gather(mleIndexZ));
    waitbar(i/round(length(zStage)/nData),hWait)
    toc
end
toc
close(hWait)
%%
if length(xStage)<length(dataOut)
    dataOut=dataOut(1:length(xStage));
    dataOutZ=dataOutZ(1:length(xStage));
end
mleBayes=(zeros(length(dataOut),3));
mleBayes(:,1)=x(ind2sub(size(x),dataOut));
mleBayes(:,2)=y(ind2sub(size(y),dataOut));
mleBayes(:,3)=z(ind2sub(size(z),dataOutZ));
% Plots
time=(1:length(dataOut))*tau;

figure(1)

plot(time,xStage(1:length(dataOut)))
hold on
plot(time,xStage(1:length(dataOut))+mleBayes(:,1))
hold off
ylabel('X Position')
legend('X Stage','X Stage plus MLE')
figure(2)
plot(time,yStage(1:length(dataOut)))
hold on
plot(time,yStage(1:length(dataOut))-mleBayes(:,2))
hold off
ylabel('Y Position')
legend('Y Stage','Y Stage plus MLE')

figure(3)
plot(time,zStage(1:length(dataOutZ)))
hold on
plot(time,zStage(1:length(dataOutZ))-mleBayes(:,3))
hold off
ylabel('Z Position')
legend('Z Stage','Z Stage plus MLE')
%%
xOut=xStage(1:length(dataOut))+mleBayes(:,1);
yOut=yStage(1:length(dataOut))-mleBayes(:,2);
zOut=zStage(1:length(dataOut))-mleBayes(:,3);

function [mleIndex,lPo]=bayesQuickCalc(lPr,lP_n_x,ktPos,n,diff,mleIndex,startIndex,nData)

for k=1:nData
    lPo=lPr+lP_n_x(:,:,ktPos(k+startIndex-1)+1,n(k+startIndex-1)+1);
    [m,mleIndex(k)]=max(lPo(:));
    lPo=lPo-m;
    %Calculate Log posterior using log-sum-exp method. This helps avoid
    %underfill. 
    lPo=lPo-log(sum(sum(exp(lPo))));
  
    
    %Account for diffusion via convolution
    lPr=log(conv2(exp(lPo),diff,'same'));
end

function [mleIndex,lPoZ]=bayesQuickCalcZ(lPrZ,lP_n_z,TAGPos,n,diffZ,mleIndex,startIndex,nData)

for k=1:nData
    lPoZ=lPrZ+lP_n_z(:,TAGPos(k+startIndex-1)+1,n(k+startIndex-1)+1)';
    [m,mleIndex(k)]=max(lPoZ);
    lPoZ=lPoZ-m;
    lPoZ=lPoZ-log(sum(exp(lPoZ))); 
    
    %Loop through, posterior becomes prior
    %Account for diffusion via convolution
    lPrZ=log(conv(exp(lPoZ),diffZ,'same'));
end

function lP_n_x=genLL(x,y,xkt,ykt,s,b,w,tau,n)
%Generate likelihood which can be indexed by position in KT and also by
%number of photons
gamma=gpuArray(zeros(size(x,1),size(x,2),length(xkt)));
lP_n_x=gpuArray(zeros(size(x,1),size(x,2),length(xkt),max(n)+1));
for i=1:length(xkt)

        %Average count rate per time      
        gamma(:,:,i)=(s*exp(-(x-xkt(i)).^2./2./w./w).*exp(-(y-ykt(i)).^2./2./w./w)+b)*tau;
        
end

for j=0:max(n)
    lP_n_x(:,:,:,j+1)=j*log(gamma)-gamma;
end

function lP_n_z=genLLZ(z,TAGPos,s,b,wZ,tau,n)
%Generate likelihood which can be indexed by position in KT and also by
%number of photons
bZ=(-0.045*z+1)*b;
bZ=b;
gamma=gpuArray(zeros(size(z,2),round(max(TAGPos))+1));
lP_n_z=gpuArray(zeros(size(z,2),round(max(TAGPos))+1,max(n)+1));
for i=0:round(max(TAGPos))

        %Average count rate per time      
        gamma(:,i+1)=(s*exp(-(z-sin(i/100+.83)).^2./2./wZ./wZ)+bZ)*tau;
        
end

for j=0:max(n)
    lP_n_z(:,:,j+1)=j*log(gamma)-gamma;
end