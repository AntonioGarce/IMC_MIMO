% This programme designs an observer zhat(t) for the state
%  in the plant used in
% SIMULINK  file newss.m. The output are the polynomials
%  nu(s),du(s), ny(s) and dy(s), where
%
%
% Zhat(s)=U(s) nu(s)/du(s)+Y(s) ny(s)/dy(s)
%
% The input is a vector with the four observer poles
%
function [nu,du,ny,dy]=lambor(P)
%
A=[0 1 0 0
-8 -6 -16 0
0 0 0 0
1 0 0 -1];
B=[0;16;0;0];C=[0 0 0 1];
J=place(A',C',P);J=J';Ao=A-J*C;
[ntu,dtu]=ss2tf(Ao,B,eye(4,4),0*B);[nty,dty]=ss2tf(Ao,J,eye(4,4),0*J);
ny=nty(1,:);ny=ny(2:5);dy=dty;du=dtu;nu=ntu(1,:);nu=nu(2:5);

