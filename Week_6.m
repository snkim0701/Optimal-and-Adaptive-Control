
%% Week_5 transfer function
% known parameters
% open loop T.R.


%% Week_5 transfer function
% unknown parameters
% open loop T.R.

clear all; clc
syms s a b k
A =[ 0 1; -a -b];
B =[0 1]';
%C = [ 1 0];
%C = [ 1 0];
C = [ 0 1];
I = eye(2,2);
TR_open = C*inv(s*I-A)*B

% close loop
FK = [k 0];
TR_closed = C*inv(s*I-(A-B*FK))*B

%% Week_5 
a12 = 1;
a22= -1;
a34= 1;
a41 = -11.65;
a43 = -a41;
b12=1;
A = [ 0 a12 0 0; 0 a22 0 0; 0 0 0 a34; a41 0 a43 0];
B = [0 b12 0 0]';
eig(A)


%% transfer function symboli math
clear all;clc
syms a b s
A =[ 0 1; -a -b];
B = [0 ; 1];
u = [ 0 1];
C= [ 1  0];
I = eye(2,2);
Kim=C*inv(s*I - (A - B*u))*B;
simplify(Kim)
%% page 49
clear all;clc
syms F M g L s k
A = [ 0 1 0 0; 0 -F/M 0 0; 0 0 0 1; -g/L 0 g/L 0];
B = [0; 1/M; 0; 0]
K= [-1 0 1 0];
C= [ -1 0 1 0];
I = eye(4,4);
Kim = (C*inv(s*I - (A-B*K))*B)  % transfer function 
                                % the ordert of tranfer function should be "4" 
simplify(Kim)
%% example 3.1 page 195
clear all; clc
syms F M g L s k1 k2 k3 k4
A = [ 0 1 0 0; 0 -F/M 0 0; 0 0 0 1; -g/L 0 g/L 0];
B = [0; 1/M; 0; 0];
K = [ k1 k2 k3 k4];
I = eye(4,4);
Den=det(s*I - (A-B*K)) % the charateristic eqn
pretty(Den)
%%
DenM = subs(Den,M,1);
DenF = subs(DenM,F,1);
DenL = subs(DenM,L,0.842);
gs = 11.65*0.842;
Deng = subs(DenL,g,gs);
pretty(Deng);

%%  .... I did not get good method to find (3-13). Try to find!

clear all; clc
syms F M g L s 
A = [ 0 1 0 0; 0 -F/M 0 0; 0 0 0 1; -g/L 0 g/L 0];
B = [0; 1/M; 0; 0];
%K = [ k1 k2 k3 k4];
K =[65.65 11.00 -72.60,-21.27]; % textbook page 196.(3-13)
I = eye(4,4);
Den=det(s*I - (A-B*K)) % the charateristic eqn
DenM = subs(Den,M,1);
DenF = subs(DenM,F,1);
DenL = subs(DenF,L,0.842);
gs = 11.65*0.842;
Deng = subs(DenL,g,gs)

%% to get (3-13) 
clear all;clc
syms s k1 k2 k3 k4
a12 = 1;
a22= -1;
a34= 1;
a41 = -11.65;
a43 = -a41;
b12=1;
A = [ 0 a12 0 0; 0 a22 0 0; 0 0 0 a34; a41 0 a43 0];
B = [0 b12 0 0]';
eig(A)
I = eye(4,4);
K = [k1 k2 k3 k4];
%K =-[65.65 11.00 -72.60,-21.27]; % textbook page 196.(3-13)
det(s*I-(A+B*K))
%
k2 =-11;
k1 = -54-233/20
k4=(108-233/20*k2+233/20)*20/233
k3=(81-233/20*k1)*20/233
K=[k1 k2 k3 k4]

% check the result
simplify(det(s*I-(A+B*K)))
eig(A+B*K)  
display('Bing go !!')
