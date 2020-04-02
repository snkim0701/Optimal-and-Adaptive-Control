

%% linear homegeneous - differential equation-ode
clear all; clc; clf
tspan = [ 0 5]; % domain
x0 = 1;         % initial point
[t,x] = ode45(@(t,x) -x, tspan,x0); % ode45 function.
plot(t,x); grid on


%% linear in-homogeneous
clear all; clc; clf
tspan = [ 0 5]; % domian
x0 = 1;         % initial point
[t,x] = ode45(@(t,x) -x+2, tspan,x0); % ode45 function.
plot(t,x); grid on

%% non-linear - converge
clear all; clc; clf
tspan = [ 0 5]; % domian
x0 = 2;         % initial point
[t,x] = ode45(@(t,x) -x - x^2 +1, tspan,x0); % ode45 function.
plot(t,x); grid on


%% non-linear - diverge
clear all; clc; clf
tspan = [ 0 0.6]; % domian
x0 = 2;         % initial point
[t,x] = ode45(@(t,x) -x + x^2 +1, tspan,x0); % ode45 function.
plot(t,x); grid on

%% (3-142) Fig 3.7  Riccati equation

clear all; clc; clf
tspan = [ 0 1]; % domian
s0 = 0;         % initial point
alpha = 0.5;
k = 150;
r = [100 1000 10000];
for i = 1:3
[t s] = ode45(@(t,s) 1 - k^2/r(i) * s^2 - 2*alpha*s, tspan,s0); % ode45 function
plot(1-t,s,'Linewidth',2); grid on; hold on
end

s0 = [0 0.19 0.5];         % initial point
r = 1000;
for i = 1:3
[t s] = ode45(@(t,s) 1 - k^2/r * s^2 - 2*alpha*s, tspan,s0(i)); % ode45 function
plot(1-t,s,'Linewidth',2); grid on; 
end
hold off
title('Fig.3.7 p(t) for various values of rho and phi')
disp('I am Fine')


%% Example 3.7 : scalar algebraic Riccati equation
syms  p a rho k

S =solve(-2*a*p+1-k^2/rho *p^2 ==0,p);
pretty(S)
disp('Bing go!')

%% Figure 3.8
clear all; clc;clf
a = 4.6;
k = 0.787;
for i=1:1000
    rho = i/10000;
    y1(i) = -(1/2)*(sqrt(a^2 + 2*k/sqrt(rho)) + sqrt(a^2 -2*k/sqrt(rho)));
    plot(real(y1), imag(y1),'Linewidth',2); hold on; grid on
    y2(i) = -(1/2)*(sqrt(a^2 + 2*k/sqrt(rho)) - sqrt(a^2 -2*k/sqrt(rho)));
    plot(real(y2), imag(y2),'Linewidth',2); 
end

%
mrho = [0.0055 0.0002 ];
for i = 1:2
    M_1(i) = -(1/2)*(sqrt(a^2 + 2*k/sqrt(mrho(i))) + sqrt(a^2 -2*k/sqrt(mrho(i))));
    plot(real(M_1), imag(M_1),'o','Markersize',10); hold on
    M_2(i) = -(1/2)*(sqrt(a^2 + 2*k/sqrt(mrho(i))) - sqrt(a^2 -2*k/sqrt(mrho(i))));
    plot(real(M_2), imag(M_2),'o', 'Markersize',10);
end
title('Loci of closed loop roots as a function of rho')



%% Figure 3.9 state -space - known parameter values

clear all; clc;clf

% First open loop system model
a = 4.6;
k = 0.787;
rho = 0.0002;
A =[0 1; 0 -a];
B = [0 k]';
C = [ 1 0];
Q =C'*C;
R = rho;

% find the optimal controller
sys = ss(A,B,C,[ ]);
[K,P,E] = lqr(sys,Q,R)

% The closed loop system
AC =A -B*B'*(1/rho)*P;
sysClosed =ss(AC,B,C,[ ]);
x0 =[0.1; 0];

[y,T,x] =initial(sysClosed,x0); grid on
subplot(2,1,1)
plot(T,y); grid on
title('the output angle with rho = 0.0002,with initial point')

u = -B'*P*(1/(rho)*x');
subplot(2,1,2)
plot(T,u); grid on
title('the control input of DC motor voltage')


%%


