function [K,J,numc,denc]=css(a,b,c,po,pc)
%
% This function compute a one d.o.f. controller for
% a n-th order SISO  STRICTLY  proper plant (continuous or discrete)
% described in state space form by the 4-tuple (a,b,c,0)
%
% 
% The user must supply a vector po with n desired
% observer poles and a  vector pc with n desired
% control poles
%
% The function returns the gains K and J and 
% the numerator and  denominator polynomial of the controller
% transfer function (see Chapter 8 in 
% "Principles of Control System Design"
% by Goodwin, Graebe and Salgado
%
%
% Warning: This function uses MATLAB function "place", thus the
%          user is advised to become familiar with the numerical
%          limitations of the algorithm used in "place".
%
[n,n]=size(a);
%
[n1,n2]=size(po);
if max(n1,n2)~=n
  error ('illegal number of observer poles')
end
 [n3,n4]=size(pc);
if max(n3,n4)~=n
  error ('illegal number of controller poles')
end
%
% Observer and controller gain computation
%
K=place(a,b,pc); J=place(a',c',po);J=J';
%
A=real(a-J*c-b*K); B=real(J);C=real(K);
%
% Controller reduction
%
[Am,Bm,Cm,Dm]=minreal(A,B,C,0);
[num,denc]=ss2tf(Am,Bm,Cm,Dm);
[nx,numc]=p_elcero(num);
%



 