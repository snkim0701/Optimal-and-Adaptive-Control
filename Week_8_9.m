
%% Week_8: https://www.mathworks.com/help/control/ref/obsv.html
clear all;clc; 
A =[1 1;4 2]; 
C = [1 0];
Q =obsv(A,C)
rank(Q)

%%
clear all;clc; 
syms a b
A =[0 1; a b]; 
C = [1 1];
Q = [C; C*A]
rank(Q)


%% observibility

A = [0 1 0;0 0 1; -1 -2 -3];
M  = eye(3);
for i = 1:3
    C = M(:,i)'
    rank(obsv(A,C))       
end 


%% The optimal observer desing: 
% Example 4.3 without observer
clear all; clc;clf
V1= 0.1;V2 = 0.01;  % noise intensity 
x0 = 1;
v1 = V1*wgn(1000,1,0);
v2 = V2*wgn(1000,1,0);
t = linspace(0,5,1000);

figure(1)
subplot(1,2,1)
plot(t,v1,'r', t,v2,'b'); grid on;
title('disturbace w1 and noise w2')
%
% system dynamics
A= -1; B = 1 ; C=1;D = 0; 
sysx = ss(A,B,C,D); 
[y,t,x]=lsim(sysx, v1, t,1);
yout = y + v2;

subplot(1,2,2)
plot(t,yout,t,y); grid on
title('the original state with disturbance and output with noise')

%the optimal observer design

% solving Riccati equation
syms Q 
eqn = 2*A*Q + V1 - Q^2/V2;
Ri = double(solve(eqn,Q))  % for Q>= 0, Ri(2) is validated.

% the optimal observer gain
ObK = Ri(2)/V2
AO= A-ObK;
BO=[1 ObK];
uO = [v1 y];
eig(AO)
sysO = ss(AO,BO,C,[ ]);
Temp=lsim(sysO,uO,t,0);
figure(3)
plot(t,yout,'b', t,Temp,'r',t,y,'k'); grid on

%% Example 5.3 - observer  design

clear all;clc;
a = 4.6;
ka = 0.787;
rho = 0.00002;
gam = 0.1;
Vd = [ 0 0; 0 gam^2]*10;
Vm = 10^-7;
A = [0 1;0 -a];
B = [ 0; ka];
C = [ 1 0];

% solve riccati equation
[X,K,L] = care(A',C',Vd,Vm)

%% 

clear all; clc
sys  = tf([10], [1 10])

bode(sys); grid on
title('bode plot of the low pass filter')

sys = ss(tf(100, [1 1 100]));

% Design LGR
K = lqry(sys,10,1) % in the cost function

% Separate control input u and the disturbance
P = sys(:,[1 1]);

% Design Kalman state estimator Kest.
Kest = kalman(P,1,0.01)

% Form LQG regulator = LQ gain + Kalman filter.
F = lqgreg(Kest,K)
% Close loop
clsys = feedback(sys,F,+1)
% Note positive feedback.

% Create the lowpass filter and add it in series with clsys.
s = tf('s');
lpf= 10/(s+10) ;
clsys_fin = lpf*clsys;

% Open- vs. closed-loop impulse responses
impulse(sys,'r--',clsys_fin,'b-')