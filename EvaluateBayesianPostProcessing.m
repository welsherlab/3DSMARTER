function EvaluateBayesianPostProcessing
%Define particle/trajectory parameters:
D=4e-12;
s=1e6; 
b=1e6;
N=5000; %approx 3 hrs comp time. 

%Define Integration time
tau=20e-6;
subtau=1e-6;
dx=1e-9; %grid spacing for Bayesian post processing position estimates.

% Use "generatetraj" to generate a 2D trajectory using Kalman tracking.

% Inputs:
% D, diffusive speed of particle in meters^2/second
% s, maximum observable count rate of particle under uniform illumination in counts/second
% b, background intensity in counts/second
% N, number of 20 us bins in original trajectory. 
% tau, bin time used for tracking in seconds. Particle position is estimated and
% used to update laser and stage positions once per tau. 
% subtau, bin time used to generate particle positions and photons in seconds

% Outputs:
% len, index of position where particle escaped tracking if applicable. 
%Otherwise len=N*number of subbins/bin.
% time, time axis of data generation in s
% stg_pos_err, error in stage position vs particle position. square root of ((x_stg-x_particle)^2+(y_stg-y_particle)^2)
% photons, number of photons generated during each subbin.
% x_particle, particle position along x for each subbin in meters.
% y_particle, particle position along y for each subbin in meters.
% x_stg, simulated stage position along x for each subbin in meters. Value only updated every tau.
% y_stage, simulated stage position along y for each subbin in meters. Value only updated each tau. 
% x_laser, laser position along x relative to the position given by x_stg in meters. 
% y_laser, laser position along y relative to the position given by y_stg in meters. 
[len,time,stg_pos_err,photons, x_part, y_part,x_stg,y_stg,x_laser,y_laser]= generatetraj(D,s,b,N,tau,subtau);
disp('done with traj generation')

% Use "analyzetrajeachsubbin" to post process the trajectory generated using generate traj in the previous step.

%Inputs:
%x_part, actual particle position along x in meters
%y_part, actual particle position along y in meters
%x_stg, piezo stage position along x in meters
%y_stg, piezo stage position along y in meters
%x_laser, laser scan center along x relative to x_stg in meters
%y_laser, laser scan center along y relative to y_stg in meters
%photons, photon counts at each subbin from trajectory
%D, assumed diffusion coefficient of particle in meters^2/second used in post processing algorithms
%subtau, time resolution of data generation in seconds used to post process.
%s, assumed count rate of particle under uniform illumination in counts/s
%b, assumed intensity of background in counts/s
%dx, spacing used to generate grid for Bayesian post processing.

%Outputs:
%Kalman_est_err, error between Kalman post processed position estimate and
%particle positions
%bayes1d_est_err, error between 1D Bayesian post processing and particle
%positions
%bayes2D_est_err, error between 2D Bayesian post processing and particle
%postions

tic
[kalman_est_err,bayes1d_est_err,bayes2d_est_err]=analyzetrajeachsubbin(x_part(1:len), y_part(1:len), x_stg(1:len), y_stg(1:len), x_laser(1:len), y_laser(1:len), photons(1:len), D, subtau,s,b,dx);
disp('done with per subbin analysis')
toc

stgpos_err=sqrt((x_part-x_stg).^2+(y_part-y_stg).^2); 
%CDF Plot
 figure
 hold on
[h,stats_stg] = cdfplot(stgpos_err);
[h,stats_k] = cdfplot(kalman_est_err);
[h,stats_1d] = cdfplot(bayes1d_est_err);
[h,stats_2d] = cdfplot(bayes2d_est_err);
hold off
legend('Stage Position','Kalman Post-Processing','1D Bayes Post-Processing','2D Bayes Post-Processing')
xlabel('Error ( m)')
ylabel('C.D.F')
ax = gca;
ax.XAxis.Exponent = -6;

%Histograms
%Manually change order of histograms so color coding still corresponds to
%CDF Plots
figure
hold on
histogram(stgpos_err,'Normalization','pdf')
histogram(kalman_est_err,'Normalization','pdf')
histogram(bayes1d_est_err,'Normalization','pdf')
histogram(bayes2d_est_err,'Normalization','pdf')
hold off
legend('Stage Position','Kalman Post-Processing','1D Bayes Post-Processing','2D Bayes Post-Processing')
xlabel('Error ( m)')
ylabel('P.D.F.')
ax = gca;
ax.XAxis.Exponent = -6;

%Particle Traj/True Postions Panel
figure
hold on
 plot(time,x_part,'Color',[19 159 255]./255,'LineWidth',2)
 plot(time,x_stg,'Color',[7 67 107]./255,'LineWidth',3)
 plot(time,y_part,'Color',[0.85 0.33 0.10],'LineWidth',2)
 plot(time,y_stg,'Color',[0.64 0.08 0.18],'LineWidth',3)
legend('x particle','x stage','y particle','y stage')
hold off
xlabel ('Time (s)')
ylabel('Position ( m)')
ax = gca;
ax.YAxis.Exponent = -6;

%Subbin photons histogram
figure
histogram(photons)
xlabel('Photons per subbin')

disp(['Every Subbin Means: kal = ' num2str(mean(kalman_est_err)) ' 1d = ' num2str(mean(bayes1d_est_err)) ' 2d = ' num2str(mean(bayes2d_est_err)) ' mean stg err = ' num2str(mean(stg_pos_err))])
disp('done')
end

function[kalman_est_err,bayes1d_est_err,bayes2d_est_err]=analyzetrajeachsubbin(x_part, y_part, x_stg, y_stg, x_laser,y_laser, photons,D,subtau,s,b,dx)

%initialize
%Generate XY coordinates
%dx=1e-10;
coord=-.25e-6:dx:.25e-6;
%2D coordinates
[x2D,y2D]=meshgrid(coord);
%1D coordinate
x1D=coord;
y1D=coord;
sigma=200e-9;
%Initialize Kalman
x_k_1=[0;0];
sigma_k_1=sigma;

%Calculate the Diff convolution matrix if needed
if D~=0
    diffull=exp(-coord.^2/4/D/subtau);
    c=find(diffull>=.1*max(diffull));
    %check symmetry of coordinates
    mid=ceil(length(diffull)/2);
    disp(['diffusion kernel is ' num2str(length(c)) ' points wide'])
    if mid-c(1)==c(end)-mid
        diffkern=diffull(c);
    else
        warning('Diffusion kernel is asymmetric')
        diffkern=diffull(c);
    end
    if length(c)<=1
        warning('diffusion matrix is delta function')
    end
    clear c mid
end

%Initialize prior
Pr2D=ones(size(x2D));       %2D
Pr1DX=ones(size(x1D));        %1D
Pr1DY=ones(size(y1D));        %1D

%Normalize Prior
Pr2D=Pr2D/sum(sum(Pr2D));       %2D
Pr1DX=Pr1DX/sum(sum(Pr1DX));    %1D
Pr1DY=Pr1DY/sum(sum(Pr1DY));    %1D

%Log Prior
lPr2D=log(Pr2D);            %2D
lPr1DX=log(Pr1DX);          %1D
lPr1DY=log(Pr1DY);          %1D

%time into signal/background estimation
testimate=15e-3;
binsestimate=int32(testimate/subtau);

for k=1:length(x_part)
    %estimate signal and bg levels based on observed photons (assume 50% s
    %and b each)
%     if k<=binsestimate
%         %use first 15 ms of data
%         guess=0.5*sum(photons(1:binsestimate))/testimate;
%         [calc_sb]=[guess;guess];
%     else
%         %use most recent 15 ms of data
%         guess=0.5*sum(photons(k-binsestimate+1:k))/testimate;
%         [calc_sb]=[guess;guess];
%     end
guess=1e6;
calc_sb=[s;b];
    
    %if signal estimate is negative or 0 crash
    if calc_sb(1)<=0
        warning('negative or 0 estimated signal')
        break
    end
    
    %if bg estimate is negative replace it with 0
    if calc_sb(2)<0
        calc_sb(2)=0;
    end
    
    %Construct likelihood of detection n photons.
    %2D
    gamma2D=calc_sb(1)*exp(-((x2D-x_laser(k)).^2./2./sigma./sigma)-((y2D-y_laser(k)).^2./2./sigma./sigma))+calc_sb(2);
    %1D
    gammaX=calc_sb(1)*exp(-(x1D-x_laser(k)).^2./2./sigma./sigma)+calc_sb(2);
    gammaY=calc_sb(1)*exp(-(y1D-y_laser(k)).^2./2./sigma./sigma)+calc_sb(2);
    
    log_n_factorial=stirling(photons(k));
    
    %2D
    lP_n_2D=photons(k)*log(gamma2D*subtau)-gamma2D*subtau-log_n_factorial;
    %1D
    lP_n_x=photons(k)*log(gammaX*subtau)-gammaX*subtau-log_n_factorial;
    lP_n_y=photons(k)*log(gammaY*subtau)-gammaY*subtau-log_n_factorial;
    
    %Calculate log posterior
    lPo2D=lPr2D+lP_n_2D; %2D
    lPoX=lPr1DX+lP_n_x; %1D
    lPoY=lPr1DY+lP_n_y; %1D
    
    %Now normalize log posterior using log-sum-exp method
    %2D
    lPo2D=lPo2D-max(max(lPo2D));
    lPo2D=lPo2D-log(sum(sum(exp(lPo2D))));
    %1D
    lPoX=lPoX-max(lPoX);
    lPoX=lPoX-log(sum(exp(lPoX)));
    
    lPoY=lPoY-max(lPoY);
    lPoY=lPoY-log(sum(exp(lPoY)));
    
    %Account for diffusion via convolution. If D=0, then posterior becomes
    %prior
    if D~=0
        % 2D conv
        lPr2D=log(conv2(diffkern,diffkern,exp(lPo2D),'same'));
        lPr2D(isinf(lPr2D))=-745;
        % 1D conv
        lPr1DX=log(conv(exp(lPoX),diffull,'same'));
        lPr1DX(isinf(lPr1DX))=-745;
        
        lPr1DY=log(conv(exp(lPoY),diffull,'same'));
        lPr1DY(isinf(lPr1DY))=-745;
    else
        lPr2D=lPo2D; %2D
        lPr1DX=lPoX; %1D
        lPr1DY=lPoY; %1D
    end
    
    %Kalman comparison
    x_k = (sigma^2*x_k_1+photons(k)*sigma_k_1*[x_laser(k);y_laser(k)])/(sigma^2+photons(k)*sigma_k_1);
    sigma_k = sigma^2*sigma_k_1/(sigma^2+photons(k)*sigma_k_1);
    
    %Update Kalman variables
    x_k_1=x_k;
    sigma_k_1=sigma_k+2*D*subtau;
    
    kalm_x(k)=x_k(1);
    kalm_y(k)=x_k(2);
    
    %Find maximum likelihood estimate of position
    %2D
    temp=x2D(lPo2D==max(max(lPo2D)));
    x_mle_2D(k)=mean(temp);
    temp=y2D(lPo2D==max(max(lPo2D)));
    y_mle_2D(k)=mean(temp);
    %1D
    temp=x1D(lPoX==max(lPoX));
    x_mle_1D(k)=mean(temp);
    temp=y1D(lPoY==max(lPoY));
    y_mle_1D(k)=mean(temp);
    
    cal_sb(k)=guess;
end

%calculate actual errors
bayes2d_est_err=sqrt((x_mle_2D+x_stg-x_part).^2+(y_mle_2D+y_stg-y_part).^2);
bayes1d_est_err=sqrt((x_mle_1D+x_stg-x_part).^2+(y_mle_1D+y_stg-y_part).^2);
kalman_est_err=sqrt((kalm_x+x_stg-x_part).^2+(kalm_y+y_stg-y_part).^2);

%make time axis so can do plotting
time=subtau:subtau:length(x_part)*subtau;
% %Make investagatory plots
% figure
% hold on
%  plot(time,x_part,'Color',[0.85 0.33 0.10],'LineWidth',2)
%  plot(time,x_stg,'Color',[0.64 0.08 0.18],'LineWidth',3)
% % plot(time,kalm_x+x_stg,'Color',[0.47 0.67 1])
% % plot(time,x_mle_2D+x_stg,'Color',[0.99,0.77,0.26])
% % plot(time,x_mle_1D+x_stg,'Color',[0.055 0.7 0])
% legend( 'Particle','Stage','Kalman','2D Bayes','1D Bayes')
% xlabel('Time (seconds)')
% title('X axis subbin estimation')
% hold off
% % 
% % figure
% % hold on
% %  plot(time,y_part,'Color',[0.85 0.33 0.10],'LineWidth',2)
% %  plot(time,y_stg,'Color',[0.64 0.08 0.18],'LineWidth',3)
% % plot(time,kalm_y+y_stg,'Color',[0.47 0.67 1])
% % plot(time,y_mle_2D+y_stg,'Color',[0.99,0.77,0.26])
% % plot(time,y_mle_1D+y_stg,'Color',[0.055 0.7 0])
% % hold off
% % legend( 'Particle','Stage','Kalman','2D Bayes','1D Bayes')
% % xlabel('Time (seconds)')
% % title('Y axis subbin estimation')
% figure
% histogram(photons)
% xlabel('photons per subbin')
% 
% figure
% histogram2([x_part-x_stg,y_part-y_stg],[kalm_x,kalm_y],'FaceColor','flat')
% view(2)
% xlabel('Actual Offset')
% ylabel('Offset Estimate')
% title('Kalman')
% 
% figure
% histogram2([x_part-x_stg,y_part-y_stg],[x_mle_1D,y_mle_1D],'FaceColor','flat')
% view(2)
% xlabel('Actual Offset')
% ylabel('Offset Estimate')
% title('1-D Bayesian')
% 
% figure
% histogram2([x_part-x_stg,y_part-y_stg],[x_mle_2D,y_mle_2D],'FaceColor','flat')
% view(2)
% xlabel('Actual Offset')
% ylabel('Offset Estimate')
% title('2-D Bayesian')
% 
% figure
% plot(time,cal_sb)
% yline(sum(photons)/(length(x_part)*subtau)/2);
% xlabel('Time (seconds)')
% ylabel('Signal/BG estimate (cps)')
% title('Signal BG estimates using subbin estimation')
% 
% figure
% hold on
% histogram(kalm_y)
% histogram(kalm_x)
% hold off
% title('per subbin histogram of raw kalman estimates')
% 
% figure
% histogram(kalman_est_err)
% title('per subbin kalman estimate error')
% 
% figure
% hold on
% histogram(y_mle_1D)
% histogram(x_mle_1D)
% title('per subbin 1D bayes MLE estimates')
% 
% figure
% histogram(bayes1d_est_err)
% title('per subbin 1D bayes estimate error')
% 
% figure
% hold on
% histogram(y_mle_2D)
% histogram(x_mle_2D)
% title('per subbin 2D bayes MLE estimates')
% 
% figure
% histogram(bayes2d_est_err)
% title('per subbin 2D bayes estimate error')
% 
% figure
% scatter(x_part-x_stg,x_mle_1D)
% xlabel('particle offset from stage')
% ylabel('position estimate')
% title('x 1d bayes')
% axis square
% xlim([-2.5e-7,2.5e-7])
% ylim([-2.5e-7,2.5e-7])
% 
% figure
% scatter(x_part-x_stg,x_mle_2D)
% xlabel('particle offset from stage')
% ylabel('position estimate')
% title('x 2d bayes')
% axis square
% xlim([-2.5e-7,2.5e-7])
% ylim([-2.5e-7,2.5e-7])
% 
% figure
% scatter(x_part-x_stg,kalm_x)
% xlabel('particle offset from stage')
% ylabel('position estimate')
% title('x kalman')
% axis square
% xlim([-2.5e-7,2.5e-7])
% ylim([-2.5e-7,2.5e-7])
% 
% figure
% scatter(y_part-y_stg,y_mle_1D)
% xlabel('particle offset from stage')
% ylabel('position estimate')
% title('y 1d bayes')
% axis square
% xlim([-2.5e-7,2.5e-7])
% ylim([-2.5e-7,2.5e-7])
% 
% figure
% scatter(y_part-y_stg,y_mle_2D)
% xlabel('particle offset from stage')
% ylabel('position estimate')
% title('y 2d bayes')
% axis square
% xlim([-2.5e-7,2.5e-7])
% ylim([-2.5e-7,2.5e-7])
% 
% figure
% scatter(y_part-y_stg,kalm_y)
% xlabel('particle offset from stage')
% ylabel('position estimate')
% title('y kalman')
% axis square
% xlim([-2.5e-7,2.5e-7])
% ylim([-2.5e-7,2.5e-7])
% 
% stepfit=table('Size',[4,3],'VariableTypes',{'double','double','double'},'RowNames',{'x slope','y slope','x y intercept', 'y y intercept'},'VariableNames',{'Kalman','1D Bayes','2D Bayes'});
% %fit x axis datas and store appropriately
% tfit([1,3],1)=polyfit(x_part-x_stg,kalm_x,1);
% tfit([1,3],2)=polyfit(x_part-x_stg,x_mle_1D,1);
% tfit([1,3],3)=polyfit(x_part-x_stg,x_mle_2D,1);
% %fit y axis datas and store appropriately
% tfit([2,4],1)=polyfit(y_part-y_stg,kalm_y,1);
% tfit([2,4],2)=polyfit(y_part-y_stg,y_mle_1D,1);
% tfit([2,4],3)=polyfit(y_part-y_stg,y_mle_2D,1);
% %disp fitting output/comparable "step test" data for position estimates
% stepfit(:,:)=array2table(tfit)
% 
% figure
% hold on
% plot(time,x_part-x_stg,'LineWidth',2)
% plot(time,kalm_x)
% plot(time,x_mle_1D,'Color',[0.4660, 0.6740, 0.1880])
% plot(time,x_mle_2D,'LineWidth',1)
% hold off
% legend('particle offset','kalman','1D Bayes','2D Bayes')
% xlabel ('time (sec)')
% ylabel ('position (m)')
% title ('x estimate comprehensive')
% 
% figure
% hold on
% plot(time,y_part-y_stg,'LineWidth',2)
% plot(time,kalm_y)
% plot(time,y_mle_1D,'Color',[0.4660, 0.6740, 0.1880])
% plot(time,y_mle_2D,'LineWidth',1)
% hold off
% legend('particle offset','kalman','1D Bayes','2D Bayes')
% xlabel ('time (sec)')
% ylabel ('position (m)')
% title ('y estimate comprehensive')

save(['scatteringproofofestimatevalue_gridsize_' num2str(dx) '.MAT'])
end

%generate trajectory using kalman position estimation and a given large bin
%time ~20e-6 usec
function [lensub,subtau_time,stg_pos_err, photonsfinegrain, x_part, y_part,stg_xout,stg_yout,laser_xout,laser_yout]=generatetraj(D,s,b,N,tau,subtau)
ki=0.032;
kp=0;

%initialize bin sizes
ratio=int32(tau/subtau);
errint_x=0;
errint_y=0;
tau_time=tau:tau:N*tau;
subtau_time=subtau:subtau:N*tau;

%generate particle positions on functionally continuous timescale
%Generate Trajectory
%Generate steps in x (dx) and y (dy)
dx = sqrt(2*D*subtau) * randn(N*ratio,1)';
dy = sqrt(2*D*subtau) * randn(N*ratio,1)';

%Generate XY trajectory
x_part = cumsum(dx);
y_part = cumsum(dy);

%PSF size
sigma=200e-9;
%w=sigma^2*eye(2,2);

%Initialize Kalman
x_k_1=[0;0];
sigma_k_1=sigma;

%Size of KT
d=500e-9;

%Stage position
prev_stg=[0;0];

% Load impulse response data
measuredImpulseResponse = dlmread('Impulse_Response.txt');
nelIR=length(measuredImpulseResponse);

%Now simulate tracking
for k=1:N
    %First, get simulated observed number of photons based on the particle
    %and laser positions
    [xlaser,ylaser]=knightsTour(k,d);   %laser positions
    laserPos=[xlaser;ylaser];
    
    %calculate observed number of photons using subtau bin sizes
    lambda=subtau*emRate(([x_part(k*ratio-ratio+1:k*ratio);y_part(k*ratio-ratio+1:k*ratio)]-prev_stg), laserPos,b,s,sigma);
    photonsfinegrain(k*ratio-ratio+1:k*ratio)=poissrnd(lambda); %number of photons for each subtau
    n=sum(photonsfinegrain(k*ratio-ratio+1:k*ratio));
    
    %Kalman comparison
    x_k = (sigma^2*x_k_1+n*sigma_k_1*laserPos)/(sigma^2+n*sigma_k_1);
    sigma_k = sigma^2*sigma_k_1/(sigma^2+n*sigma_k_1);
    
    %Update Kalman variables
    x_k_1=x_k;
    sigma_k_1=sigma_k+2*D*tau;
    
    kalm_x(k)=x_k(1);
    kalm_y(k)=x_k(2);
    
    %Integrate errors and Move "stage"
    %integrate error
    errint_x=errint_x+x_k(1);
    errint_y=errint_y+x_k(2);
    
    %Apply integrated error times integral weighting constant to recenter
    %particle
    stageComX(k)=kp*kalm_x(k)+ki*errint_x;
    stageComY(k)=kp*kalm_y(k)+ki*errint_y;
    if k<=nelIR
        stg_x=impulseResponse(stageComX,measuredImpulseResponse);
        stg_y=impulseResponse(stageComY,measuredImpulseResponse);
    else
        stg_x=impulseResponse(stageComX(end-nelIR+1:end),measuredImpulseResponse);
        stg_y=impulseResponse(stageComY(end-nelIR+1:end),measuredImpulseResponse);
    end
    stg_xout(k*ratio-ratio+1:k*ratio)=stg_x(end);
    stg_yout(k*ratio-ratio+1:k*ratio)=stg_y(end);
    prev_stg=[stg_x(end);stg_y(end)];
    
    laser_xout(k*ratio-ratio+1:k*ratio)=xlaser;
    laser_yout(k*ratio-ratio+1:k*ratio)=ylaser;
    
    %if particle escaped box end trajectory
    if abs(prev_stg(1)-x_part(k*ratio))>.5e-6
        disp('particle escaped box along x')
        break
    end
    
    if abs(prev_stg(2)-y_part(k*ratio))>.5e-6
        disp('particle escaped box along y')
        break
    end
end
%% Calculate Remaining Outputs and Plot
len=k;
lensub=k*ratio;
% sanitize err calc inputs
asdf=size(stg_xout);
if asdf(1)~=1
    stg_xout=stg_xout';
    stg_yout=stg_yout';
end

% calculate err
stg_pos_err=sqrt((stg_xout(1:lensub)-x_part(1:lensub)).^2+(stg_yout(1:lensub)-y_part(1:lensub)).^2);

% %% Make Plots
% % plot particle positions, stage positions, and postition estimates
% % relative to stage position on same plot. 
% figure()
% clf
% hold on
% plot(subtau_time,x_part,'Color',[0.07,0.62,1],'LineWidth',2)
% plot(subtau_time(1:len*ratio),stg_xout,'Color',[0 0 1],'LineWidth',3)
% plot(tau_time(1:len),kalm_x,'Color',[0.47 0.67 1])
% plot(subtau_time,y_part,'Color',[0.85 0.33 0.10],'LineWidth',2)
% plot(subtau_time(1:len*ratio),stg_yout,'Color',[0.64 0.08 0.18],'LineWidth',3)
% plot(tau_time(1:len),kalm_y,'Color',[0.99,0.77,0.26])
% hold off
% legend( 'particle x','stage x','kalman x','particle y','stage y','kalman y')
% xlabel('Time (seconds)')

% %Plot histogram of kalman position estmates relative to center of stage
% %position.
% figure
% hold on
% histogram(kalm_x)
% histogram(kalm_y)
% hold off
% title('kalman position estimates from tau during initial traj generation')

% plot kalman position estimates on trajectory generation tau. 
% kalmesterr=sqrt((stg_xout(1:lensub)+kalm_x-x_part).^2+(stg_yout(1:lensub)+kalm_y-y_part).^2);
% figure
% histogram(kalmesterr)
% title('err in kalman estimate during traj gen')

end

%Stirling approximation for log(n!)
function out=stirling(n)
if n<10
    out=log(factorial(n));
else
    out=n.*log(n)-n+1;
end
end

%4pt scan pattern coordinates
function [xOut,yOut]=knightsTour(index,d)
x=d*[0.5;-0.5;-0.5;0.5];
y=d*[-0.5;-0.5;0.5;0.5];
xOut=x(mod(index-1,length(x))+1);
yOut=y(mod(index-1,length(x))+1);
end

% Impulse Response Function
function out = impulseResponse(x,measuredImpulseResponse)
%Convolve with stage command data
temp = conv(x,measuredImpulseResponse);
%Crop data
out = temp(1:length(x));
end

% Photon Count Rate at a Given Stage Position
function out=emRate(x,c,b,s,sigma)
%x - position of particle
%c - position of beam
%b - background
%w - spatial covariance is calculated from sigma (1/sigma^2)*eye
%sigma - laser standard deviation
%code now uses expanded/parallelized version of below which is equivalent
%to original implementation to 15 decimal places and 63x faster (for subbin:bin
%ratio of 2000)
%out=s*exp(-0.5*(x-c)'*inv(w)*(x-c))+b;
%out=s*exp(-0.5*(x-c)'/w*(x-c))+b;
out=s*exp(((-0.5/sigma^2)*(x(1,:)-c(1)).^2)+((-0.5/sigma^2)*(x(2,:)-c(2)).^2))+b;
end
