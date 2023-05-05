 function [X,Y] = paq(A,B,C) 
%
% FUNCTION [X,Y] = paq(A,B,C)
%
% Program to solve a pole assignment equation using
% least squares fitting in the frequency domain.
% this program is suitable for continuous time and for discrete
% time in delta form ( it is not numerically 
% suitable for discrete time in z form)
%
%       A*X+B*Y=C
%
% The user must provide polynomials A and C, in MONIC FORM,
% and polynomial B.
%
%     DEGREE OF B MUST BE LESS THAN DEGREE
%     OF A, AND DEGREE OF A MUST BE LESS OR  EQUAL THAN DEGREE OF C.
%
%     A,B and C must be coprime polynomials
%
% The program returns MONIC polynomial X and
%  (not necessarily  monic) polynomial Y
%
%
%           The  directory where the 
%           function p_elcero.m resides
%           must be  in the working path

%
% Verification of monic character of A and C
%

%
%First eliminate leading zeros
%
[ntrue,A]=p_elcero(A);
[ntrue,B]=p_elcero(B);
[ntrue,C]=p_elcero(C);
%
%
if A(1)~=1
	error('A must be a monic polynomial')
end
if C(1)~=1
	error('C must be a monic polynomial')
end

% Degrees and verification
la=length(A);lb=length(B);lc=length(C);
da=la-1;db=lb-1;dc=lc-1;dx=dc-da;lx=dx+1;
%
%number of coefficients (excluding the monic coeff) 
na=da;nb=lb;nc=dc;nx=dx;ny=na;dy=ny-1;
if db >= da
	error('degree of B  must be less than degree of A')
end
if da > dc
	error('degree of A  must be less or equal than degree of C')
end
%
%Checking (quasi) cancellations 
%
[n1,d1]=minreal(A,C,.001);
   [ntrue,n1]=p_elcero(n1);
   if length(n1)<length(A)
     error (' (Quasi) cancellation of factors between A and C. Simplify first')
   end
[n1,d1]=minreal(B,C,.001);
   [ntrue,n1]=p_elcero(n1);
   if length(n1)<length(B)
     error (' (Quasi) cancellation of factors between B and C. Simplify first')
   end
[n1,d1]=minreal(B,A,.001);
   [ntrue,n1]=p_elcero(n1);
   if length(n1)<length(B)
     error (' (Quasi) cancellation of factors between A and B. Simplify first')
   end
%
%
%
% Building the frequency vector
%
%
 wmax=2*max ( (C(lc))^(1/nc), C(2)/nc);
%
dw=wmax/1000;dwl=log10(dw);w=logspace(dwl,log10(wmax),1000);
%
% Initial vectors and matrices
%
if nx>0
fx=ones(nx,1);
end
fy=ones(ny,1); P=10000*eye(nx+ny);theta=zeros(nx+ny,1);
%
j=sqrt(-1);
%
% Recursive least squares computation of the coeffs of X and Y
%
for i=1:1000
 v=j*w(i);

%
 Mc=polyval(C,v)-polyval(A,v)*v^dx;
   if nx>0 
     for k=dx-1:-1:0
       fx(dx-k)=v^k;
     end
   end
   for k=dy:-1:0
    fy(dy-k+1)=v^k;
   end
 if nx>0
   fi=[fx*polyval(A,v);fy*polyval(B,v)];
  else
   fi=fy*polyval(B,v);
 end
fir=real(fi);K=P*fir;D=1+fir'*P*fir;M=real(Mc);
theta=theta+K*(M-fir'*theta)/D;
P=P-P*fir*fir'*P/D;
fii=imag(fi);K=P*fii;D=1+fii'*P*fii;M=imag(Mc);
theta=theta+K*(M-fii'*theta)/D;
P=P-P*fii*fii'*P/D;
end
%
%
if nx>0
 X=[1;theta(1:nx,1)];
 else
 X=1;
end
Y=theta(nx+1:nx+ny,1);








