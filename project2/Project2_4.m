% EC 3310 Computer Assignment #2 Part4
% Alessio Grompone 02/03/2013
% Target Motion with 1 sensors and Weighted Least Squares Algorithm

%DATA
nloops=100;

delta=0.5;    %sample time
nsamples=5; %number of samples

S1p=[0;0];      %Sensor 1 position
sigma1r=100;   %standard deviation sensor 1 range
sigma1b=1*pi/180;   %standard deviation sensor 1 bearing
S1v=diag([sigma1r^2;sigma1b^2]);  %covariance
S1=[sigma1r;sigma1b];


ht=135*pi/180; %Target Attitude (from degrees to radians) 
st=360*2000/3600;   %Target Speed (from knots to yards/sec)
%Inital Conditions
xo=10000;
yo=10000;
Vxo=st*sin(ht);
Vyo=st*cos(ht);
%State Vector
xi=[xo;Vxo;yo;Vyo];

%Transition Matrix
F=[1,delta,0,0;0,1,0,0;0,0,1,delta;0,0,0,1];

%Measurement Matrix
H=[1,0,0,0;
    0,0,1,0];

for kk=1:nloops,
%initialize Matixes
zout=[];   %measurements output
posout=[]; %true target position
error=[];  %distance errors
polar=[];
xe=[];
S1vc=[];
x=xi;
y=[];


for ii=1:nsamples,
    
    ztrue=H*x;
    posout=[posout,ztrue];
    
    %% SENSOR 1
   
    range=sqrt((ztrue-S1p)'*(ztrue-S1p));
    bearing=atan2(ztrue(2)-S1p(2),ztrue(1)-S1p(1));
    polar=[range;bearing];    
    rrr=randn(size(S1));
    zpol=polar+rrr.*S1;%has measurements adding random error
    xx=zpol(1)*cos(zpol(2));
    yy=zpol(1)*sin(zpol(2));
    cartesian=[xx;yy]+S1p;  %cartesian coordinates
    
    %collects the real position values in one matrix
    zpolar=[zpolar,zpol];
    zout=[zout,cartesian];  %collects the measurement values in one matrix
    ztilde=ztrue-(cartesian); %error between real position and measured position
    error=[error,sqrt(ztilde'*ztilde)];%collects error values
   
    x=F*x;
 
%Estimating covarinace in cartesian
   Fx=[cos(zpol(2)), -zpol(1)*sin(zpol(2));
       sin(zpol(2)), zpol(1)*cos(zpol(2))];
   
       S1vcv=Fx*S1v*Fx';
       
       y=[y;cartesian];

end
 %State vector estimate
    %(Forward Projection Method)
   
    K=F;
    M=[H;H*K;H*K^2;H*K^3;H*K^4];

   W=inv([S1vcv,zeros(2),zeros(2),zeros(2),zeros(2);
   zeros(2),S1vcv,zeros(2),zeros(2),zeros(2);
   zeros(2),zeros(2),S1vcv,zeros(2),zeros(2);
   zeros(2),zeros(2),zeros(2),S1vcv,zeros(2);
   zeros(2),zeros(2),zeros(2),zeros(2),S1vcv]);
    
    
    xe1=inv(M'*W*M)*M'*W*y; 
    xe2=K*xe1;
    xe3=K*K*xe1;
    xe4=K*K*K*xe1;
    xe5=K*K*K*K*xe1;
    
    Xe=[xe1,xe2,xe3,xe4,xe5];
    
    Pe=H*Xe; 
Pmerror=[];

Perror=posout-Pe; 
%Pmerror=sqrt(diag((Perror'*Perror)));
for zi=1:nsamples,
    
Pmerror=[Pmerror,sqrt(Perror(:,zi)'*Perror(:,zi))];

end   
    
    
if kk==1
    merror=error;
    zoutmean=zout;
    xmean=Pe;
    Pmean=Pmerror;
    x45=xe5;
    
else
    merror=merror+error;
    zoutmean=zoutmean+zout; 
    xmean=xmean+Pe;
    Pmean=Pmean+Pmerror;
    x45=x45+xe5;
   
end
end

merror=merror/nloops;
zoutmean=zoutmean/nloops;
xmean=xmean/nloops;
Pmean=Pmean/nloops;
x45=x45/nloops

%plot1 Target motion with mean noise measurements
plot(posout(1,:),posout(2,:),'-',zoutmean(1,:),zoutmean(2,:),'-');
title('Target Motion and Mean Measurement Trajectory');
% title('Target Motion and Mean Measurement Trajectory');
pause;

%plot1 Target motion with mean state vector estimations
plot(posout(1,:),posout(2,:),'-',xmean(1,:),xmean(2,:),'-');
title('Target Motion and Least Square Estimates');
% title('Target Motion and Mean Measurement Trajectory');
pause;

% %Distance errors between Target and measurements
time=[1:max(size(Pmean))]*(2*delta)-(2*delta);
plot(time,Pmean,'-');
hold
plot(time,merror,'-');
title('Distance Error in Estimation vs Time (100 runs)');
pause;
