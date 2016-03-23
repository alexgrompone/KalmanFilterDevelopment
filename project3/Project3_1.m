% EC 3310 Computer Assignment #3 Part1
% Alessio Grompone 02/03/2013
% Kalman Filter SWAG
clear all
%DATA
nloops=100;

delta=0.1;    %sample time
nsamples=1000; %number of samples

q=0; %plant noise

S1p=[0;0];      %Sensor 1 position
sigma1r=100;   %standard deviation sensor 1 range
sigma1b=1*pi/180;   %standard deviation sensor 1 bearing
S1v=diag([sigma1r^2;sigma1b^2]);  %covariance
S1=[sigma1r;sigma1b];

%Inital Conditions

xo=-6000;
yo=20000;
Vxo=100;
Vyo=-173.2051;
%State Vector
xi=[xo;Vxo;yo;Vyo];

%Transition Matrix
F=[1,delta,0,0;0,1,0,0;0,0,1,delta;0,0,0,1];

%Measurement Matrix
H=[1,0,0,0;
    0,0,1,0];
%SWAG initialization
xos=0;
yos=10000;
Vxos=0;
Vyos=-500*2/3.6;
xo_k=[xos;Vxos;yos;Vyos];
Po_k=(10^10)*eye(4);

%Kalman Filter Parameters

Qk=(q)*[(delta^3)/3,(delta^2)/2,0,0;
    (delta^2)/2,delta,0,0;
    0,0,(delta^3)/3,(delta^2)/2;
    0,0,(delta^2)/2,delta];

for kk=1:nloops,
%initialize Matixes
zout=[];   %measurements output
posout=[]; %true target position
error=[];  %distance errors
polar=[];
x=xi;
Xk_k=xo_k;
Pk_k=Po_k;
Pest=[];
Perror=[];
MaxEigs=[];
Rk11=[];
for ii=1:nsamples,
    
    ztrue=H*x;
    posout=[posout,ztrue];
    
    %% SENSOR 1
   
    range=sqrt((ztrue-S1p)'*(ztrue-S1p));
    bearing=atan2(ztrue(2)-S1p(2),ztrue(1)-S1p(1));
    polar=[range;bearing];    
    rrr=randn(size(S1));
    z=polar+rrr.*S1;%has measurements adding random error
    xpol=z(1)*cos(z(2));
    ypol=z(1)*sin(z(2));
    cartesian=[xpol;ypol]+S1p;  %cartesian coordinates
    
    %collects the real position values in one matrix
    zout=[zout,cartesian];  %collects the measurement values in one matrix
    ztilde=ztrue-(cartesian); %error between real position and measured position
    error=[error,sqrt(ztilde'*ztilde)];%collects error values
   
    
%Estimating covarinace in cartesian
Fx=[cos(z(2)), -z(1)*sin(z(2));...
       sin(z(2)), z(1)*cos(z(2))];  

Rk1=Fx*S1v*Fx';
Rk11=[Rk11,Rk1];
  
%measurement Ztilda an K
if ii ==1
    Xk1_k = xo_k;
    Pk1_k = Po_k;
end
Ztk1=cartesian-H*Xk1_k; 
Kk1=Pk1_k*H'*inv(H*Pk1_k*H' + Rk1); 

%correction
Xk1_k1=Xk1_k + Kk1*(Ztk1);
Pk1_k1=(eye(4)-Kk1*H)*Pk1_k*(eye(4)-Kk1*H)' + Kk1*Rk1*Kk1';

%Reset estimate & cov for the next cycle
Xk_k = Xk1_k1;
Pk_k = Pk1_k1;

%prediction
x=F*x;
Xk1_k=F*Xk_k;
Pk1_k=F*Pk_k*F' + Qk;

%Position Estimates and Errors from the filter
Pestout=H*Xk_k;
Pest=[Pest,Pestout];
Ptilde=ztrue-Pestout; %error between real position and measured position
Perror=[Perror,sqrt(Ptilde'*Ptilde)];%collects error values

MaxEig=sqrt(eigs(Pk_k,1));

if ii<11
MaxEigs=[MaxEig];
Perror2=[sqrt(Ptilde'*Ptilde)];
else
MaxEigs=[MaxEigs;MaxEig]; 
Perror2=[Perror2,sqrt(Ptilde'*Ptilde)];
end

end

if kk==1
    merror=error;
    zoutmean=zout;
    mXest=Pest;
    Pmerror=Perror;
    mMaxEigs=MaxEigs;
    Pmerror2=Perror2;
else
    merror=merror+error;
    zoutmean=zoutmean+zout;  
    mXest=mXest+Pest;
    Pmerror=Pmerror+Perror;
    mMaxEigs=mMaxEigs+MaxEigs;
    Pmerror2=Pmerror2+Perror2;
end
end

merror=merror/nloops;
zoutmean=zoutmean/nloops;
mXest=mXest/nloops;
Pmerror=Pmerror/nloops;
mMaxEigs=mMaxEigs/nloops;
Pmerror2=Pmerror2/nloops;

%plot1 Target motion with mean noise measurements
plot(posout(1,:),posout(2,:),'-',zoutmean(1,:),zoutmean(2,:),'-');
title('Target Motion and Mean Measurement Trajectory');
% title('Target Motion and Mean Measurement Trajectory');
pause;

%plot2 Target motion with Mean Track Trajectory
plot(posout(1,:),posout(2,:),'-',mXest(1,:),mXest(2,:),'-');
title('Target Motion and Mean Track Trajectory');
pause;

%Distance errors between Mean Track Trajectory and measurements
time=[1:max(size(error))]*delta-delta;
plot(time,merror,'-');
hold
time=[1:max(size(Pmerror))]*delta-delta;
plot(time+0.1,Pmerror,'-','Color','red');
title('Mean Measurements Errors and Mean track position errors vs Time (100 runs)');
pause;
hold off

%Distance errors between Mean Track Trajectory and Max Eigenvalues
time=[1:max(size(MaxEigs))]*delta-delta;
plot(time,mMaxEigs,'-');
hold
plot(time,Pmerror2,'-','Color','red');
title('MaximumEigenvalues Errors vs Mean Tracking Errors(100 runs)');
pause;
hold off
close all