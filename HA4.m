

clear all; clc;clf

% system dynamics
a= 4.6;
gam = 0.1;
Vdg =gam^2*10;            % the intensity of the disturbance
Vm =10^-7;                % the noise intensity 
%ka = 0.787;
% the plant state space model

A = [ 0 1;0 -a];
B = [0 1]';
C = [1 0];
sys= ss(A,B,C,0); 

% an trajectory example with an initial point
x0 =[1;1];
y0 = C*x0;  % inital values of the plant states 
N = 1000;  % the sampling number of data
tf  = 10 ;    % the final time 
t = linspace(0,tf,N);
vd = Vdg*randn(N,1);   % gaussian noise ... if awgn is available, use agwn (it stands for 
                       %   additive white gaussian noise   
[y,t,x]=lsim(sys,zeros(N,1), t,x0);   % without disturnbance 
[yd,t,x]=lsim(sys,vd, t,x0);          % disturnbance 

figure(1)
plot(t,y,'r', t,yd,'b'); grid on
title('without disturbace  and with disturbance')
axis([0 tf 0.8 2])

% the plant output corrupted by a noise
vm = Vm*randn(N,1);
ym = yd + vm;
figure(2)
plot(t,yd,'b', t,ym, 'r'); grid on
axis([0 tf 0.8 2])
title('without(blue) and with (red) measurement noise')
%%  design the optimal observer gain
[L,P,E] = lqe(A,B,C,Vdg,Vm)
%%
% Total system 
AT  = [A  zeros(2,2) ; L*C  A-L*C];
BT =  [ B  B zeros(2,2); zeros(2,2) B L];
CT = [C zeros(1,2); zeros(2,2) eye(2,2)];

sysT = ss(AT,BT,CT,0);
uT = [zeros(N,1),vd ,zeros(N,1), vm];

[yTd,t,x] = lsim(sysT,  uT, t,[x0; 0;0]); % input with  u = 0 and  disturnbance 
yTdd = yTd(:,1) + vm;
figure(3)
plot(t,yTdd,'k',t,yTd(:,2),'b', t,yTd(:,3),'r'  ); grid on ;
title('observe the x1 (black and blue) and x2(with red)')
axis([0 10  -10   20])

%%
% figure(3)
% plot(diff(yTd(:,1)),'b'); hold on
% plot(yTd(:,3),'r')