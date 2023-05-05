clc;
clear;
close all;

syms s x1 x2 x3 u u1 u2 u3 x_punkt x1_punkt x2_punkt x3_punkt;
%% Modellparameter
a1 = 0.00751;
a2 = 0.00418;
a3 = 0.05;
a4 = 0.03755;
a5 = 0.02091;
a6 = 0.00315;

b1 = 0.00192;
b2 = 0.05;
b3 = 0.00959;
b4 = 0.1866;
b5 = 0.14;

k1 = 0.01061;
k2 = 2.5;
k3 = 6.84;
k4 = 2.5531;

SP = [a1, a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,k1,k2,k3,k4];
%% nonlinear System
x1_punkt = a1*x3 + a2*x2 - b1*u1 -b2*u2 - k1;
x2_punkt = -a3*x2*u2 + k2;
x3_punkt = -a4*x3 - a5*x2 +b3*u1 + ((a6*x3+b4)/(b5*u3+k3))*u3 + k4;
y1 = x1;
y2 = x2;
y3 = x3;

%% workingpoint
xR = [1; 15; 70];
uR = [214.13; 3.33; 65.40-0.02];
x0 = [1; 25; 50];
uR1 = [214.13; 3.33+0.02; 65.40];

%% Linearization
x = [x1 x2 x3];
y = [y1 y2 y3];
u = [u1 u2 u3];
n = 3;                                      %Variablen
m = 3;                                      %Inputs
p = 3;                                      %Outputs
x_punkt = [x1_punkt x2_punkt x3_punkt]'; 
linear_A = jacobian(x_punkt, x);
linear_B = jacobian(x_punkt, u);
linear_C = jacobian(y, x);

A= subs(linear_A,x',xR);
A= subs(A,u',uR);
A = double(A);

B= subs(linear_B,x',xR);
B= subs(B,u',uR);
B = double(B);

C= subs(linear_C,x',xR);
C= subs(C,u',uR);
C = double(C);

D = zeros(p,m);
sys = ss(A,B,C,D);
%% Stabilität des AP
eig(A);
%% Steuerbarkeit des AP (controllability of working point)
Qs = [B A*B];
RangQs = rank(Qs);
%% Beobachtbarkeit des AP (observability)
Qb = [C; C*A]
RangQb = rank(Qb)
%% StateSpaceSystem
ZRM = ss(A, B, C, D ) 
%% Null- & Polstellen
zero (ZRM)
pole (ZRM)
%% Übertragungsfunkion (transferfunction)
G_tf = minreal (tf ( ZRM ))
%% Übertragungsfunktion Symbolisch (transferfunction symbolic)
G_sym = simplify (C*inv (s*eye (3, 3) - A)*B + D);
pretty ( G_sym )
%% Inversion
G_star_sym = inv ( G_sym );
pretty ( G_star_sym )
%% Kontrolle (controll)
ident = simplify ( G_sym * G_star_sym )
%% Propering Filter
T_filt_sym = 0.5e-1;
G_filt_sym = 1/( T_filt_sym *s + 1);
G_star_filt_sym = G_star_sym * G_filt_sym ;
pretty ( G_star_filt_sym )
%% Übertagungsfunktion (transferfunction)
G_star_filt_tf = sym2tf ( G_star_filt_sym )
%% ZRM (statespacemodell)
G_star_filt_ss = ss ( G_star_filt_tf )
[a,b,c,d] = ssdata(G_star_filt_ss)
%% Verstärker (gain)
k = 1e11;

K = [k 0, 0;
     0, k, 0;
     0, 0, k];

xR_punkt = originalModel(xR,uR);
uR_punkt_offset = inv(B)*xR_punkt;
uR_offset = uR - uR_punkt_offset;
%% Funktion definieren
function G_tf = sym2tf ( G_sym )
% SYM2TF Symbolic transfer function matrix to numerical transfer function matrix .
%
%   SYM2TF ( G_SYM ) returns the normalized numerical transfer function matrix
%   representation of the symbolic transfer function matrix G_SYM .
%
%   Example :
%
%       sym2tf ([s/(s+1) , (s +2)/(2* s +1)])
%
%       returns
%
%       Transfer function from input 1 to output :
%       s
%       -----
%       s + 1
%
%       Transfer function from input 2 to output :
%       0.5 s + 1
%       ---------
%       s + 0.5
%
% Copyright Joerg J. Buchholz , Hochschule Bremen , buchholz@hs - bremen .de
% Determine the numbers of rows and columns of the symbolic transfer function matrix
[ n_rows , n_cols ] = size ( G_sym );
% Disassemble every single symbolic transfer function into numerator and denominator
[ num_sym , den_sym ] = numden ( G_sym );
% Loop over all rows
for i_row = 1 : n_rows
    % Loop over all columns
    for i_col = 1 : n_cols
        % Transform the symbolic numerator of the current transfer function
        % to numerical ( coefficients of the polynomial )
        num_tf {i_row , i_col } = sym2poly ( num_sym (i_row , i_col ));
        % Transform the symbolic denominator of the current transfer function
        % to numerical ( coefficients of the polynomial )
        den_tf {i_row , i_col } = sym2poly ( den_sym (i_row , i_col ));
        % Normalize , so that leading denominator coefficient equals 1
        num_tf {i_row , i_col } = num_tf {i_row , i_col }/ den_tf {i_row , i_col }(1);
        den_tf {i_row , i_col } = den_tf {i_row , i_col }/ den_tf {i_row , i_col }(1);
    end
end
% Assemble the numerical transfer function matrix
G_tf = tf ( num_tf , den_tf );
end

function x_punkt = originalModel(x, u)
    a1 = 0.00751;
    a2 = 0.00418;
    a3 = 0.05;
    a4 = 0.03755;
    a5 = 0.02091;
    a6 = 0.00315;
    
    b1 = 0.00192;
    b2 = 0.05;
    b3 = 0.00959;
    b4 = 0.1866;
    b5 = 0.14;
    
    k1 = 0.01061;
    k2 = 2.5;
    k3 = 6.84;
    k4 = 2.5531;

    x1_punkt = a1*x(3)+a2*x(2)-b1*u(1)-b2*u(2)-k1;
    x2_punkt = -a3*x(2)*u(2)+k2;
    x3_punkt = -a4*x(3)-a5*x(2)+b3*u(1)-(a6*x(3)+b4)*u(3)/(b5*u(3)+k3)+k4;
    
    x_punkt = [x1_punkt;x2_punkt;x3_punkt];
end


    