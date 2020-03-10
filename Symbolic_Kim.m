
%% https://www.mathworks.com/help/symbolic/performing-symbolic-computations.html
% differentiation
% one variable

syms x
f = sin(x)^2;
diff(f)
%% two variables
syms x y
f = sin(x)^2 + cos(y)^2;
diff(f,x)
diff(f,y)

%% Solve Algebraic Equations 

% One symbolic variables
syms x
solve(x+3 ==0)
solve(x^2 +3*x +2 ==0)

%% multivariables
clear all; clc;clf
syms f(x,y)
f(x,y) =2*x^2 + 2*x*y + 4*y^2;
fcontour(f,'Levellist',[1 4 10],'Linewidth',2); grid on
axis([-3 3 -3 3]);
hold on

g(x,y) = x+ y;
fcontour(g,'LevelList',[2])


%% symbolic math - linear equations
syms x y
f = x + 2*y;
g =2*x +y;
eqn =[f == 5, g == 4];
S =solve(eqn,x,y);
S.x
S.y
%%  undefined linear equations
syms x y a
f = a*x + 2*y;
g =2*x +y;
eqn =[f == 5, g == 4];
S =solve(eqn,x,y);
S.x
S.y
%% parameter optimization - Bryson
clear all; clc;
syms x u k a b m c % z  = lambda 
f = x+m*u-c ;
H = (1/2) *(x^2/a^2 + u^2/b^2) + k*f;
eqns =[f==0, diff(H,x)==0,diff(H,u) ==0];

S = solve(eqns,x,u,k);
sol=[S.x; S.u; S.k]  
%% optimal problems 
clear all; clc;
syms x u k a b m c % z  = lambda 
f = x+m*u-c ;
eqns =[f==0,x/a^2 + k ==0,u/b^2 + k*m ==0];
S = solve(eqns,x,u,k);
[S.x;S.u,;S.k]




