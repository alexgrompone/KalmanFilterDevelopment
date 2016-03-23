% EC 3310 Computer Assignment #5 Part1
% Alessio Grompone 02/03/2013
% Extended Kalman Filter
clear all
%DATA
nloops=100;
delta=0.1;    %sample time
nsamples=600; %number of samples
q2=4;         %plant noise square

S1p=[0;0];      %Sensor 1 position
sigma1r=100;    %standard deviation sensor 1 range
sigma1b=1*pi/180;   %standard deviation sensor 1 bearing
S1v=[sigma1r^2 0;0 sigma1b^2];  %covariance
S1=[sigma1r;sigma1b];

%Inital Conditions
xo=6000;
Vxo=100;
yo=8000;
Vyo=-173.2051;
%State Vector
xi=[xo;Vxo;yo;Vyo];
%Transition Matrix
F=[1,delta,0,0;0,1,0,0;0,0,1,delta;0,0,0,1];
%Measurement Matrix
H=[1,0,0,0;
    0,0,1,0];

Qko=(q2)*[(delta^3)/3,(delta^2)/2,0,0;
    (delta^2)/2,delta,0,0;
    0,0,(delta^3)/3,(delta^2)/2;
    0,0,(delta^2)/2,delta];

for kk=1:nloops, 
%initialize Matrixes
zout=[];   %measurements output
posout=[]; %true target position
error=[];  %distance errors
polar=[];
x=xi;
Pest=[];
Perror=[];
MaxEigs=[];
ChikV=[];
q2V=[];
q2VK=[];
home=[];
Rate=[];
homeerr=[];
for ii=1:nsamples,
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %%%%%%%%%%%%%%%%%% True Trajectory %%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qk=(q2)*[(delta^3)/3,(delta^2)/2,0,0,0;
    (delta^2)/2,delta,0,0,0;
    0,0,(delta^3)/3,(delta^2)/2,0;
    0,0,(delta^2)/2,delta,0;
   
    0,0,0,0,(0.01)*delta];   % 
  % 0,0,0,0,(0.001)*delta];   % a  
  % 0,0,0,0,((0.0001))*delta]; % b

    a=-60; %yards/sec^2
    v=sqrt(x(2)^2+x(4)^2);
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
    cartesian=[xpol;ypol]+S1p;%cartesian coordinates
    
    %collects the real position values in one matrix
    zout=[zout,cartesian];  %collects the measurement values in one matrix
    ztilde=abs(ztrue)-abs(cartesian); %error between real position and measured position
    error=[error,sqrt(ztilde'*ztilde)];%collects error values
    
    Fx=[cos(z(2)), -z(1)*sin(z(2));
        sin(z(2)),  z(1)*cos(z(2))];  
    Rk1=Fx*S1v*Fx';%Estimating covarinace in cartesian 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Extended Kalman FILTER %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=[1 0 0 0;1/delta 0 -1/delta 0; 0 1 0 0; 0 1/delta 0 -1/delta];  
if ii==1 %momorization of first two measurements and R
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
if ii<3 %stop filter for the first two measurements
Pest=[];
Ptilde=0;
Perror=[];
MaxEig=0;
Xk1_k=[1 0 0 0; 1/delta 0 -1/delta 0; 0 1 0 0; 0 1/delta 0 -1/delta]*[z2;z1];
Pk1_k=10^10*eye(4);
Xk1_k(5)=0.01;
Pk1_k(5,5)=0.001;
else
if ii==3
Xk1_k=[1 0 0 0; 1/delta 0 -1/delta 0; 0 1 0 0; 0 1/delta 0 -1/delta]*[z2;z1];
Pk1_k=M*[Rk12 zeros(2,2); zeros(2,2) Rk11+H*inv(F)*Qko*(inv(F))'*H']*M';
Xk1_k(5)=0.01;
Pk1_k(5,5)=0.001;
end

Xxf=[Xk1_k(1)+((sin(Xk1_k(5)*delta))*Xk1_k(2))/Xk1_k(5)-((1-cos(Xk1_k(5)*delta))*Xk1_k(4))/Xk1_k(5);
        cos(Xk1_k(5)*delta)*Xk1_k(2)-sin(Xk1_k(5)*delta)*Xk1_k(4);
        ((1-cos(Xk1_k(5)*delta))*Xk1_k(2))/Xk1_k(5)+Xk1_k(3)+((sin(Xk1_k(5)*delta))*Xk1_k(4))/Xk1_k(5);
        sin(Xk1_k(5)*delta)*Xk1_k(2)+cos(Xk1_k(5)*delta)*Xk1_k(4);
        Xk1_k(5)];

fo1=(((Xk1_k(4)*(1-cos(Xk1_k(5)*delta)))+Xk1_k(5)*delta*(((cos(Xk1_k(5)*delta))*Xk1_k(2))-((sin(Xk1_k(5)*delta))*Xk1_k(4))))/(Xk1_k(5)^2))-(((sin(Xk1_k(5)*delta))*Xk1_k(2))/(Xk1_k(5)^2));
fo2=-delta*(Xk1_k(2)*sin(Xk1_k(5)*delta)+Xk1_k(4)*cos(Xk1_k(5)*delta));
fo3=(((Xk1_k(2)*(cos(Xk1_k(5)*delta)-1))+Xk1_k(5)*delta*(((cos(Xk1_k(5)*delta))*Xk1_k(4))+((sin(Xk1_k(5)*delta))*Xk1_k(2))))/(Xk1_k(5)^2))-(((sin(Xk1_k(5)*delta))*Xk1_k(4))/(Xk1_k(5)^2));
fo4=delta*(Xk1_k(2)*cos(Xk1_k(5)*delta)-Xk1_k(4)*sin(Xk1_k(5)*delta));

Ftt=[1     sin(Xk1_k(5)*delta)/Xk1_k(5)               0     -(2*(sin(Xk1_k(5)*delta/2))^2)/Xk1_k(5)    fo1;
     0     cos(Xk1_k(5)*delta)                  0     -sin(Xk1_k(5)*delta)                 fo2;
     0     (2*((sin(Xk1_k(5)*delta/2))^2))/Xk1_k(5)   1     (sin(Xk1_k(5)*delta))/Xk1_k(5)             fo3;
     0     sin(Xk1_k(5)*delta)                  0     cos(Xk1_k(5)*delta)                  fo4;
     0           0                        0            0                       1];
%Pk1_k=Ftt*Pk1_k*Ftt'+ Qk;%%%
%Measurement Update
sqx=(Xxf(1)^2)+(Xxf(3)^2);
Hk=[Xxf(1)/(sqrt(sqx)) 0 Xxf(3)/(sqrt(sqx)) 0 0;
    -Xxf(3)/sqx 0 Xxf(1)/sqx 0 0];

Kk1=Pk1_k*Hk'*inv(Hk*Pk1_k*Hk' + S1v); 
Ztk1=z-[sqrt(sqx);atan2(Xxf(3),Xxf(1))]; 
%chi square shift trigger
% Chik=Ztk1'*inv(Hk*Pk1_k*Hk' + S1v)*Ztk1;
% ChikV=[ChikV,Chik];
% 
% % chis=abs(2*pi*(Hk*Pk1_k*Hk' + S1v))^(-0.5)*expm(-0.5*(Ztk1'*inv(Hk*Pk1_k*Hk' + S1v)*Ztk1));
% if Chik>=10.6
%     q2=1000;
% elseif Chik<=1.5
%     q2=1;
% else
%     q2=q2;
% end

%Correction 
Xk1_k1=Xxf + Kk1*(Ztk1);
Pk1_k1=(eye(5)-Kk1*Hk)*Pk1_k*(eye(5)-Kk1*Hk)' + Kk1*S1v*Kk1';
Pestoutkf=[Xk1_k1(1);Xk1_k1(3)];
Pest=[Pest,Pestoutkf];
Ptildekf=ztrue-Pestoutkf; %error between real position and measured position
Perror=[Perror,sqrt(Ptildekf'*Ptildekf)];%collects error values
MaxEig=sqrt(eigs(Pk1_k1,1));
if ii<11
MaxEigs=[MaxEig];
Perror2=[sqrt(Ptildekf'*Ptildekf)];
else
MaxEigs=[MaxEigs;MaxEig]; 
Perror2=[Perror2,sqrt(Ptildekf'*Ptildekf)];
end
if ii<250
    x=F*x;
    Rat=0;
elseif ii<350
    hom=(a/v);
    A=[0 1 0 0; 0 0 0 -hom; 0 0 0 1; 0 hom 0 0];
    Ft=expm(A*delta);
    x=Ft*x;
    Rat=hom;
else
    x=F*x;
    Rat=0;
end
  %Reset estimate & cov for the next cycle
%Xk1_k = Ftt*Xk1_k1;
Xk1_k = Xk1_k1;
Pk1_k=Ftt*Pk1_k1*Ftt'+ Qk;
%Pk1_k=Pk1_k1;
home=[home,Xk1_k(5)];
homeerr=[homeerr,abs(Xk1_k(5)-Rat)];
Rate=[Rate,Rat];
end % if statement for two measurement initialization
end

if kk==1
    merror=error;
    zoutmean=zout;
    mXest=Pest;
    Pmerror=Perror;
    mMaxEigs=MaxEigs;
    Pmerror2=Perror2; 
    Omega=home;
    homeerrV=homeerr;
else
    merror=merror+error;
    zoutmean=zoutmean+zout;  
    mXest=mXest+Pest;
    Pmerror=Pmerror+Perror;
    mMaxEigs=mMaxEigs+MaxEigs;
    Pmerror2=Pmerror2+Perror2; 
    Omega=Omega+home;
    homeerrV=homeerrV+homeerr;
end
end

merror=merror/nloops;
zoutmean=zoutmean/nloops;
mXest=mXest/nloops;
Pmerror=Pmerror/nloops;
mMaxEigs=mMaxEigs/nloops;
Pmerror2=Pmerror2/nloops;
Omega=Omega/nloops;
homeerrV=homeerrV/nloops;

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
plot(time,Pmerror2,'-','Color','red');
title('Maximum Eigenvalues Errors vs Mean Tracking Errors (100 runs)');
pause;
hold off
close all

%Mean estimated turnrate and real turnrate
time=(1:max(size(Omega)))*delta-delta;
plot(time+0.2,Omega,'-');
hold
time=(1:max(size(Rate)))*delta-delta;
plot(time,Rate,'-','Color','red');
title('Mean estimated turnrate and real turnrate (100 runs)');
pause;
hold off
close all

%Mean estimated turnrate error
time=(1:max(size(homeerrV)))*delta-delta;
plot(time+0.2,homeerrV,'-');
title('Mean turn rate error (100 runs)');
pause;
hold off
close all