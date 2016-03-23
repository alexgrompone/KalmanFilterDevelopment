% EC 3310 Computer Assignment #1 (part1)
% Alessio Grompone 01/18/2013
% Target Motion

%DATA

delta=0.1;    %sample time
nsamples=250; %number of samples
sigmarange=100;   %standard deviation on range
sigmaalfa=3*pi/180;   %standard deviation on angle
sigma=[sigmarange;sigmaalfa];  %standard deviations in one vector
ht=126.8699*pi/180; %Target Attitude (from degrees to radians) 
st=360*2000/3600;   %Target Speed (from knots to yards/sec)
%Inital COnditions
xo=6000;
yo=8000;
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
%initialize Matixes
zout=[];   %measurements output
posout=[]; %true target position
error=[];  %distance errors
polar=[];
x=xi;

for ii=1:nsamples,
    ztrue=H*x;
    posout=[posout,ztrue]; %collects the real position values in one matrix
    range=sqrt(ztrue(1)^2+ztrue(2)^2);
    bearing=atan2(ztrue(2),ztrue(1));
    polar=[range;bearing];
   
z=polar+randn(size(sigma)).*sigma;%has measurements adding random error
    xpol=z(1)*cos(z(2));
    ypol=z(1)*sin(z(2));
    cartesian=[xpol;ypol];
    zout=[zout,cartesian];%collects the measurement values in one matrix
    ztilde=ztrue-cartesian; %error between real position and measured position
    error=[error,sqrt(ztilde'*ztilde)];%collects error values
    x=F*x; %target motion ( x2 = TransMatrix * x1 )
end

%plot1 Target motion and true position
 
 plot(posout(1,:),posout(2,:),'-');
 title('True Target Motion');
 axis([5500 10500 4000 9000]);
pause
%plot2 Target motion with noise measurements

plot(posout(1,:),posout(2,:),'-',zout(1,:),zout(2,:),'-');
title('Target Motion and Measurements per Sample Times');
axis([5500 10500 4000 9000]);
pause;

%Distance errors between Target and measurements
time=[1:max(size(zout))]*delta;
plot(time,error,'-');
title('Distance Error in Measurements vs Time');
axis([0 25 0 1500]);




