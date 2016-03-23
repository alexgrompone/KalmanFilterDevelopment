% EC 3310 Computer Assignment #2 Part 1
% Alessio Grompone 01/24/2013
% Target Motion with 2 sensors

%DATA
nloops=100;

delta=0.5;    %sample time
nsamples=5; %number of samples

S1p=[0;0];      %Sensor 1 position
sigma1r=100;   %standard deviation sensor 1 range
sigma1b=3*pi/180;   %standard deviation sensor 1 bearing
S1v=diag([sigma1r^2;sigma1b^2]);  %covariance
S1=[sigma1r;sigma1b];

S2p=[20000;0];  %Sensor 2 position
sigma2r=100;   %standard deviation sensor 2 range
sigma2b=3*pi/180;   %standard deviation sensor 2 bearing
S2v=diag([sigma2r^2;sigma2b^2]);  %covariance
S2=[sigma2r;sigma2b];

ht=135*pi/180; %Target Attitude (from degrees to radians) 
st=360*2000/3600;   %Target Speed (from knots to yards/sec)
%Inital COnditions
xo=10000;
yo=10000;
Vxo=st*sin(ht);
Vyo=st*cos(ht);
%State Vector
xi=[xo;Vxo;
    yo;Vyo];
%Transition Matrix
F=[1,delta,0,0;
    0,1,0,0;
    0,0,1,delta;
    0,0,0,1];
%Measurement Matrix
H=[1,0,0,0;
    0,0,1,0];

for kk=1:nloops, 
%initialize Matixes
zout=[];   %measurements output
posout=[]; %true target position
error=[];  %distance errors
polar=[];
x=xi;

for ii=1:nsamples,
    
    ztrue=H*x;
    posout=[posout,ztrue];
    
    %% SENSOR 1
   
    range=sqrt((ztrue-S1p)'*(ztrue-S1p));
    bearing=atan2(ztrue(2)-S1p(2),ztrue(1)-S1p(1));
    polar=[range;bearing];    
    rrr=randn(size(S1));
    z=polar+rrr.*S1;%has measurements adding random error
    xx=z(1)*cos(z(2));
    yy=z(1)*sin(z(2));
    cartesian=[xx;yy]+S1p;  %cartesian coordinates
    
    %collects the real position values in one matrix
    zout=[zout,cartesian];  %collects the measurement values in one matrix
    ztilde=ztrue-(cartesian); %error between real position and measured position
    error=[error,sqrt(ztilde'*ztilde)];%collects error values
    
    x=F*x;
%     %% SENSOR 2
    range=sqrt((ztrue-S2p)'*(ztrue-S2p));
    bearing=atan2(ztrue(2)-S2p(2),ztrue(1)-S2p(1));
    polar=[range;bearing];    
    rrr=randn(size(S2));
    z=polar+rrr.*S2;%has measurements adding random error
    xx=z(1)*cos(z(2));
    yy=z(1)*sin(z(2));
    cartesian=[xx;yy]+S2p; %cartesian coordinates
    
    %collects the real position values in one matrix
    zout=[zout,cartesian];  %collects the measurement values in one matrix
    ztilde=ztrue-(cartesian); %error between real position and measured position
    error=[error,sqrt(ztilde'*ztilde)];%collects error values

     x=F*x; %target motion ( x2 = TransMatrix * x1 )
    
end

if kk==1
    merror=error;
    zoutmean=zout;
else
    merror=merror+error;
    zoutmean=zoutmean+zout; 
end
end

merror=merror/nloops;
zoutmean=zoutmean/nloops;

%plot1 Target motion with noise measurements
plot(posout(1,:),posout(2,:),'-',zoutmean(1,:),zoutmean(2,:),'-');
title('Target Motion and Mean Measurement Trajectory');

pause;

%Distance errors between Target and measurements
time=[1:max(size(error))]*delta-delta;
plot(time,merror,'-');
title('Distance Error in Measurements vs Time (100 runs)');
 