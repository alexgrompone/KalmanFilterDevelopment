%EC 3310 Grompone Project #1

%TargetMotion

delta=0.1;    %sample time

nsamples=250; %number of samples


%standard deviations
sigmax=100;
sigmay=100;
sigma=[sigmax;sigmay];
%Target parameters
ht=126.8699*pi/180; %Attitude
st=360*2000/3600;   %Speed
%State Vector
xi=[6000;st*sin(ht);8000;st*cos(ht)];
%Transition Matrix
F=[1,delta,0,0;
    0,1,0,0;
    0,0,1,delta;
    0,0,0,1];
%Measurement Matrix
H=[1,0,0,0;0,0,1,0];
%initialize Matixes
zout=[]; %output
posout=[]; %true target position
error=[]; %distance errors
x=xi;

for ii=1:nsamples,
    ztrue=H*x;
    posout=[posout,ztrue];
    z=ztrue+randn(size(sigma)).*sigma;
    zout=[zout,z];
    ztilde=ztrue-z;
    error=[error,sqrt(ztilde'*ztilde)];
    x=F*x;
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



