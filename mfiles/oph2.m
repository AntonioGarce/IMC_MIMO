%function [numq,denq,ne,de]=oph2(numvv,denvv,numww,denww)
%
%  FUNCTION [numq,denq,ne,de]=oph2(numww,denww,numvv,denvv)
%
%  Function to perform H2 minimization to solve the
%  model matching problem. 
%  i.e. to find a rational stable transfer function Q_o(s)
%  such that
%
% Q_o(s)=arg min ||W(s)-V(s)Q(s)||_2
%
% where W(s)=numww/denww and V(s)=numvv/denvv are rational 
% stable tranfer functions, with W(s) possibly NMP and
%  V(s) definitely having zeros in the open RHP and, possibly,
% in the open LHP.
%
% The function returns the numerator of Q_o(s), numq, 
% and the denominator of Q_o(s), denq.
% It also returns numerator, ne, and denominator, de,
% of the error function W(s)-V(s)Q_o(s)
%
% Note that Q_o(s) results, in general to be improper
%
% For further reference see:
% Control System Design (Goodwin, Graebe and Salgado)
% and
% Feedback Control Theory  (Doyle, Francis and Tannenbaum)
%


%
% Factorise V(s) in an all pass function and a minimum phase 
% function
rvv=roots(numvv);nvv=length(rvv);
%
k=0;l=0;
for i=1:nvv
 if real(rvv(i))>0
      rvvap(k+1)=rvv(i);k=k+1; % NMP root
%    elseif real(rvv(i))==0
%      error('vv(s) has a zero on the imaginary axis')
    else
      rvvmp(l+1)=rvv(i);l=l+1; % LHP root 
 end
end
%
if k==0
  error( 'vv(s) has no NMP zeros')
end
%
  nvvmp=numvv(1)*poly(-rvvap);
if l~=0
 nvvmp=conv(poly(rvvmp), nvvmp);
end
%
dvvmp=denvv;    
%
nvvap=poly(rvvap);dvvap=poly(-rvvap);
%
nx=conv(numww,dvvap);dx=conv(denww,nvvap);
[axx,bxx,cxx,dxx]=tf2ss(nx,dx);
%
[a1,b1,c1,d1,a2,b2,c2,d2,m]=stabproj(axx,bxx,cxx,dxx);
%
[ny,dy]=ss2tf(a1,b1,c1,d1,1);
%
numqx=conv(ny,dvvmp);denqx=conv(dy,nvvmp);
numqx=numqx/denqx(1);denqx=denqx/denqx(1);
[nt,numqx]=p_elcero(numqx);
 if length(numqx)<=length(denqx)
    [numq,denq]=minreal(numqx,denqx);
 else
    [denq,numq]=minreal(denqx,numqx);[nt,denq]=p_elcero(denq);
 end
numq=numq/denq(1);denq=denq/denq(1);
%
[nep,dep]=parallel(numww,denww,-conv(numvv,numq),conv(denvv,denq));
[nt,nep]=p_elcero(nep);
[ne,de]=minreal(nep,dep,.01);


