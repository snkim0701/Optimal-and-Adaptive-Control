

%%   Figure 3.9 state -space - known parameter values

clear all; clc;clf

% First open loop system model
a = 4.6;
k = 0.787;

A =[0 1;0 -a];
B = [0 k]';
C = [ 1 0];

%  the open loop transfer function
[num,den]= ss2tf(A,B,C,0);
eig(A)
x0 =[0.1; 0];
sysO = ss(A,B,C,[ ]);


tf = 10;     % the final time
N = 1000; 
t = linspace(0, tf,N);  % the number of the sampling
u = zeros(N,1);

% initial point response
u = zeros(size(t));
[y,t,x]=lsim(sysO, u, t,x0); 

figure(1)
subplot(2,1,1)
plot(t, y,'r'); grid on
axis([-0.1  tf -1.5 2])
title('without input , the output due to the initial points')
xlabel('time')
ylabel('output')

% step response
u = zeros(size(t));
u(t>0) = 1;
[y,t,x]=lsim(sysO, u, t,x0);  

subplot(2,1,2)
plot(t,u,'b',t, x,'r'); grid on
title('step input(blue) and the states(red)')
xlabel('time')
ylabel('output')

%% the LQR controller

% find the optimal controller
rho = 0.00002;
Q =C'*C;
R = rho;
[K,P,E] = lqr(sysO,Q,R)

% The closed loop system
AC =A -B*B'*(1/rho)*P;
sysClosed =ss(AC,B,C,[ ]);
[num,den]= ss2tf(AC,B,C,0);

% eig(AC)

% initial point response
x0 =[0.1  0];
u = zeros(size(t));
[y,t,x]=lsim(sysClosed, u, t,x0);      

figure(2)
subplot(2,1,1)
plot(t, y,'r'); grid on
axis([ -0.1  0.5 -0.1 0.11]) ; grid on
title('without input , the output due to the initial points')
xlabel('time')
ylabel('output')
ux = -K*x';
subplot(2,1,2)
plot(t, ux,'r'); grid on
axis([ -0.1  0.5 -25  5]) ; grid on
title('the control input')
xlabel('time')
ylabel('input')

%  the step input response
setPoint = zeros(size(t));
setPoint(t>0) = 1;
figure(3)
% setPoint = sin(2*pi*t);
plot(t,setPoint,'b'); hold on; grid on
axis([ -0.1  tf  -1   2]) ; 
x0 = [ 0 0]';
%
Hc = 0.787/175.9785  % at s = 0. the closed loop transfer function gain
[y,t,x]=lsim(sysClosed,Hc^-1*setPoint, t,x0);     
% subplot(2,1,2)
plot(t, y,'r')
axis([ -0.1  tf  -1   2]) ; grid on
title('set point (reference point)  reference(t) = 1')
xlabel('time')
ylabel('LQR output with unit step input')

%% P controller


a = 4.6;
k = 0.787;
PG =6.7217;

AP = [0 1  ; -PG*k -a  ];
% eig(AP)

BP = [ 0  PG*k]';
C = [ 1 0 ];
sysP = ss(AP,BP,C,0)
[num den] = ss2tf(AP,BP,C,0);

%  P response
tf = 100;     % the final time
N = 1000; 
t = linspace(0, tf,N);  % the number of the sampling
u = zeros(size(t));
x0 =[0  0];
[y,t,x]=lsim(sysP, u, t,x0);      

figure(4)
subplot(2,1,1)
plot(t, y,'r'); grid on
axis([ -0.1  tf -0.1 1]) ; grid on
title('without input , the output due to the initial points')
xlabel('time')
ylabel('PI-output')

% the step input response
setPoint = zeros(size(t));
setPoint(t>0) = 1;
%
[y,t,x]=lsim(sysP, setPoint, t,x0);      
subplot(2,1,2)
plot(t, y,'r')

figure(5)
plot(t, y,'r')
axis([ -0.1  tf -1 2]) ; grid on
title('set point (reference point)  reference(t) = 1')
xlabel('time')
ylabel('P output with unit step input')





%% PI 

% the characteristic equation of the closed loop

syms PG IG s

% a = 4.6;
% k = 0.787;
% API = [0 1 0 ; -PG*k -a  k; -IG 0 0];
% I = eye(3,3);
% det(s*I -API)

% real poles 

PG = 8.96;
IG = 4.58;
API = [0 1 0 ; -PG*k -a  k; -IG 0 0];
%eig(API)
BP = [ 0  PG*k  IG]';
C = [ 1 0 0];
sysPI = ss(API,BP,C,0);
% [num den] = ss2tf(API,BP,C,0);
%  PI response

% initial responce 

% the step input response
setPoint = zeros(size(t));
setPoint(t>0) = 1;
x0 = [ 0 0 0]';
[y,t,x]=lsim(sysPI, setPoint, t,x0);      
figure(6)
plot(t, y,'r')
axis([ -0.1  tf -1 2]) ; grid on
title('set point (reference point)  reference(t) = 1')
xlabel('time')
ylabel('PI output with unit step input')




















%%

a = 4.6;
k = 0.787;
PG = 8.9581;
IG = 4.5743;
AP = [0 1 0 ; -PG*k -a  k; -IG 0 0];
% eig(AP)

BP = [ 0  PG*k  IG]';
C = [ 1 0 0];
sysPI = ss(AP,BP,C,0)
[num den] = ss2tf(AP,BP,C,0);
%% PI response

u = zeros(size(t));
x0_PI = [x0  0]';
[y,t,x]=lsim(sysPI, u, t,x0_PI);      

figure(4)
subplot(2,1,1)
plot(t, y,'r'); grid on
axis([ -0.1  tf -0.1 1]) ; grid on
title('without input , the output due to the initial points')
xlabel('time')
ylabel('PI-output')

% the step input response
setPoint = zeros(size(t));
setPoint(t>0) = 1;
%
[y,t,x]=lsim(sysPI, setPoint, t,x0_PI);      
subplot(2,1,2)
plot(t, y,'r')
axis([ -0.1  0.5 -1 2]) ; grid on
title('set point (reference point)  reference(t) = 1')
xlabel('time')
ylabel('LQR output with unit step input')




%%




plot(T,u); grid on
title('the control input of DC motor voltage')





%% 


clear all; clc
% system paramters 

a = 4.6;
ka = 0.787;
% system dynamics
A = [ 0 1;0 -a];
B = [0 k]';
C = [1 0];



gam = 0.1;
Vdg=gam^2*10;            % the intensity of the disturbance
Vm =10^-7;                     % the noise intensity 
ka = 0.787;
% the plant state space model
sys= ss(A,B,C,0); 

%% initial point response

x0 =[-1;-1];
y0 = C*x0;  % inital values of the plant states 
N = 1000;  % the sampling number of data
tf  = 50 ;    % the final time 
t = linspace(0,tf,N);
vd = Vdg*randn(N,1);  % gaussian noise ... if awgn is available, use agwn (it stands for 
                                      %   additive white gaussian noise        
[y,t,x]=lsim(sys, zeros(N,1), t,x0);   % without disturnbance                                      
[yd,t,x]=lsim(sys,  vd, t,x0);              % disturnbance 

figure(1)
plot(t,y,'r', t,yd,'b'); grid on
title('without disturbace  and with disturbance')
axis([0 tf -1.5 1])
%% step response
step = ( t>=0);
% figure(2)
% plot(t,step)
% axis([0  tf -1.5 5])
%%
x0 = [0;0];
[y,t,x]=lsim(sys, step, t,x0);   % without disturnbance                                      
[yd,t,x]=lsim(sys, step+vd, t,x0);              % disturnbance 

figure(3)
plot(t,y,'r', t,yd,'b'); grid on
title('without disturbace  and with disturbance')
axis([0 tf -0.1 2])

%% PID controller
p = 4;
Ac = A - p*B*C;
eig(Ac)
sysP = ss(Ac,B,C,0);
stepinfo(sys)
stepinfo(sysP)
%%
[y,t,x]=lsim(sysP, step, t,x0);   % without disturnbance                                      
[yd,t,x]=lsim(sysP, step+vd, t,x0);              % disturnbance 

figure(3)
plot(t,y,'r', t,yd,'b'); grid on
title('without disturbace  and with disturbance')
axis([0 tf -0.1 2])

%%
Q = eye(2,2);
R = 0.0001;
[K,S,E] = lqr(sys,Q,R);

Aclqr = A - B*K;
sysLQR  = ss(Aclqr,B,C,0);
[y,t,x]=lsim(sysLQR, step, t,x0);   % without disturnbance                                      
[yd,t,x]=lsim(sysLQR, step+vd, t,x0);              % disturnbance 

figure(4)
plot(t,y,'r', t,yd,'b'); grid on
title('without disturbace  and with disturbance')
axis([0 tf -0.1 2])
stepinfo(sysLQR)
%%   
B = [ 0 ka]';
[L,P,E] = lqe(A,B,C,Vdg,Vm)
%%

% an trajectory example with an initial point
x0 =[-1;-1];
y0 = C*x0;  % inital values of the plant states 
N = 1000;  % the sampling number of data
tf  = 5 ;    % the final time 
t = linspace(0,tf,N);
vd = Vdg*randn(N,1);  % gaussian noise ... if awgn is available, use agwn (it stands for 
                                      %   additive white gaussian noise        
[y,t,x]=lsim(sys, zeros(N,1), t,x0);   % without disturnbance                                      
[yd,t,x]=lsim(sys,  vd, t,x0);              % disturnbance 

figure(1)
plot(t,y,'r', t,yd,'b'); grid on
title('without disturbace  and with disturbance')
axis([0 tf -1.5 0])

%%  the plant output corrupted by a noise
vm = Vm*randn(N,1);
ym = yd + vm;
figure(2)
plot(t,yd,'b', t,ym, 'r'); grid on
axis([0 tf -1.5 0])
% axis([0 10 -1 1.5])
title('with disturbance and noise')

%% 
%  design the optimal observer gain
[L,P,E] = lqe(A,B,C,Vdg,Vm)


%%
% Total system 
AT  = [A  zeros(2,2) ; L*C  A-L*C];
BT =  [ B  B zeros(2,2); zeros(2,2) B L];
CT = [C zeros(1,2); zeros(2,2) eye(2,2)];

u =zeros(N,1);
% u = sin(2*pi*t);
sysT = ss(AT,BT,CT,0);
%uT = [zeros(N,1),vd ,zeros(N,1), vm];
uT = [u,vd ,u, vm];
[yTd,t,x] = lsim(sysT,  uT, t,[x0; 0;0]); % input with  u = 0 and  disturnbance 
yTdd = yTd(:,1) + vm;
figure(3)
plot(t,yTdd,'k',t,yTd(:,2),'b', t,yTd(:,3),'r'  ); grid on ; hold on
title('without input. The output(black), The x1 observer(blue), the x2 observer(red)')
axis([0 tf  -15  5])

%%  Total system 
% AT  = [A  zeros(2,2) ; L*C  A-L*C];
% BT =  [ B  B zeros(2,2); zeros(2,2) B L];
% CT = [C zeros(1,2); zeros(2,2) eye(2,2)];
% sysT = ss(AT,BT,CT,0);

u = sin(2*pi*t);
uT = [ka*u,vd ,ka*u, vm];
[yTd,t,x] = lsim(sysT,  uT, t,[x0; 0;0]); % input with  u = 0 and  disturnbance 
yTdd = yTd(:,1) + vm;
figure(3)
plot(t,yTdd,'k',t,yTd(:,2),'b', t,yTd(:,3),'r'  ); grid on ; hold on
title('with input.')
axis([0 tf  -15  5])


%%
eig(AT)
eig(A)
eig(A-L*C)
%%
TFs = tf(sysO)











