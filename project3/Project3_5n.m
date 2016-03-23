% EC 3310 Computer Assignment #3 Part5
% Alessio Grompone 02/03/2013
% Kalman Filter with Target Maneuver
clear all
%DATA
nloops=100;

delta=0.1;    %sample time
nsamples=600; %number of samples

q2=1000; %plant noise square

S1p=[0;0];      %Sensor 1 position
sigma1r=100;   %standard deviation sensor 1 range
sigma1b=1*pi/180;   %standard deviation sensor 1 bearing
S1v=diag([sigma1r^2;sigma1b^2]);  %covariance
S1=[sigma1r;sigma1b];

%Inital Conditions

xo=6000;
yo=8000;
Vxo=100;
Vyo=-173.2051;
%State Vector
xi=[xo;Vxo;yo;Vyo];

%Transition Matrix
F=[1,delta,0,0;0,1,0,0;0,0,1,delta;0,0,0,1];

%Measurement Matrix
H=[1,0,0,0;
    0,0,1,0];

%Kalman Filter Parameters

Qk=(q2)*[(delta^3)/3,(delta^2)/2,0,0;
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
Pest=[];
Perror=[];
MaxEigs=[];


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
    ztilde=abs(ztrue)-abs(cartesian); %error between real position and measured position
    error=[error,sqrt(ztilde'*ztilde)];%collects error values
   

    if ii<250
    x=F*x;
    elseif ii<350 
    a=-60; %yards/sec^2
    v=sqrt(x(2)^2+x(4)^2);
    %hom=((a/v)*180)/pi;   
    hom=(a/v); 
    A=[0 1 0 0; 0 0 0 -hom; 0 0 0 1; 0 hom 0 0];
    Ft=expm(A*delta);    
    x=Ft*x;
    else
    x=F*x;  
    end
%Estimating covarinace in cartesian  
Fx=[cos(z(2)), -z(1)*sin(z(2));
       sin(z(2)), z(1)*cos(z(2))];  

Rk1=Fx*S1v*Fx';

%measurement Initialized FILTER

M=[1 0 0 0;1/delta 0 -1/delta 0; 0 1 0 0; 0 1/delta 0 -1/delta];  

%momorization of first two measurements and R
if ii==1
    z1=cartesian;
    z2=cartesian;
    Rk11=Rk1;
    Rk12=Rk1;
    
elseif ii==2
    z1=z1;
    z2=cartesian;
    Rk11=Rk11;
    Rk12=Rk1;
else
    z1=z1;
    z2=z2;
    Rk11=Rk11;
    Rk12=Rk12;
end

%stop filter for the first two measurements
if ii<3
Pest=[];
Ptilde=0;
Perror=[];
MaxEig=0;
Xk_k=[0;0;0;0];
Pk_k=eye(4,4);

else
%Initailization
if ii==3
Xk1_k=[1 0 0 0; 1/delta 0 -1/delta 0; 0 1 0 0; 0 1/delta 0 -1/delta]*[z2;z1];

Pk1_k=M*[Rk12 zeros(2,2); zeros(2,2) Rk11+H*inv(F)*Qk*(inv(F))'*H']*M';
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
pause;

%plot2 Target motion with Mean Track Trajectory
plot(posout(1,:),posout(2,:),'-',mXest(1,:),mXest(2,:),'-');
title('Target Motion and Mean Track Trajectory');
pause;

%Distance errors between Mean Track Trajectory and measurements
time=[1:max(size(error))]*delta-delta;
plot(time,merror,'-');
hold
time=[1:max(size(Perror))]*delta-delta;
plot(time+0.2,Pmerror,'-','Color','red');
title('Mean Measurements Errors and Mean track position errors vs Time (100 runs)');
pause;
hold off

%Distance errors between Mean Track Trajectory and Max Eigenvalues
time=(1:max(size(MaxEigs)))*delta-delta;
plot(time,mMaxEigs,'-');
hold
plot(time,Pmerror2,'-','Color','black');
title('Maximum Eigenvalues Errors vs Mean Tracking Errors (100 runs)');
pause;
hold off
close all