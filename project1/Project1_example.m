% EC 3310 Computer Assignment #1 (part1)
% Alessio Grompone 01/18/2013
% Target Motion

%DATA

delta=0.1;    %sample time
nsamples=250; %number of samples
sigmax=100;   %standard deviation along x
sigmay=100;   %standard deviation along y
sigma=[sigmax;sigmay];  %standard deviations in one vector
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
x=xi;

for ii=1:nsamples,
    ztrue=H*x;
    posout=[posout,ztrue]; %collects the real position values in one matrix
    z=ztrue+randn(size(sigma)).*sigma;%has measurements adding random error
    zout=[zout,z];  %collects the measurement values in one matrix
    ztilde=ztrue-z; %error between real position and measured position
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
title('Target Motion and Measurements at Simple Times');
axis([5500 10500 4000 9000]);
pause;


%Distance errors between Target and measurements
time=[1:max(size(zout))]*delta;
plot(time,error,'-');
title('Error in Measurements vs Time');
axis([0 25 0 400]);



