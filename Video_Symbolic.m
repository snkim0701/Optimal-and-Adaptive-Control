
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

%% Home assignment_1
clear all; clc;
syms x y z la1 la2
f1 = x +2*y +3*z -10;
f2 = x -y +2*z -1;
H = x^2 + y^2 + z^2 + la1*f1 + la2*f2;
eqns =[f1==0, f2==0, diff(H,x) == 0, diff(H,y) ==0, diff(H,z)==0];
S= solve(eqns,x,y,z, la1,la2)
[S.x,S.y,S.z]

%% Home assignment_2
clear all; clc;
syms x y a b la
f = x^2/(a^2) + y^2/(b^2) -1;
H =4*(x+y) +la*f;
eqns =[f==0,diff(H,x) == 0, diff(H,y) ==0];
S= solve(eqns,x,y,la)
[S.x,S.y]
4*(S.x + S.y)
simplify(4*(S.x + S.y))

%% Home assignment_3 
clear all; clc;
syms x y z a b c la
f = x^2/(a^2) + y^2/(b^2) + z^2/(c^2)-1;
H =8*(x*y*z) +la*f;
eqns =[f==0,diff(H,x) == 0, diff(H,y) ==0,diff(H,z)==0];
S= solve(eqns,x,y,z,la);
[S.x,S.y,S.z]

%% Transfer function
% problem_1  and problem_2
clear all; clc
syms s a b k
A =[ 0 1; -a -b];
B =[0 1]';
C = [ 1 0];
%C = [ 0 1];
I = eye(2,2);
TR_open = C*inv(s*I-A)*B

% close loop
FK = [k 0];
TR_closed = C*inv(s*I-(A-B*FK))*B

%% Transfer function
% problem_3 
clear all; clc
syms s a b k1 k2
A =[ 0 1; -a -b];
B =[0 1]';
C = [ 1 0];
%C = [ 0 1];
I = eye(2,2);
TR_open = C*inv(s*I-A)*B

% close loop
FK = [k1 k2];
TR_closed = C*inv(s*I-(A-B*FK))*B

%% Transfer function
% problem_3 
clear all; clc
syms s a b k1 k2
A =[ 0 1; -a -b];
B =[0 1]';
C = [ 1 1];
%C = [ 0 1];
I = eye(2,2);
TR_open = C*inv(s*I-A)*B

% close loop
FK = [k1 k2];
TR_closed = C*inv(s*I-(A-B*FK))*B