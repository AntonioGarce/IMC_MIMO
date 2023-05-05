  function sm=smax(p,z,epsi,wl,wh)
%
%  FUNCTION sm=smax(p,z,epsi,wl,wh)
%
%  Function to compute a lower bound for the
%  peak of the  nominal  sensitity So.
%
%  This is done in the following setup:
%
%  1)The plant model has a number of unstable poles, 
%    and the effect of one particular   zero in the
%    open RHP is examined. 
%
%  2)Design specifications require that
%
%    |So|<epsi  for omega<wl
%    |To|<epsi  for omega >wh>wl
%
%  The user must provide the following data:
%
%  p : vector with  RHP poles (p=0 if no unstable poles exists)
%  z : NMP  zero
%  epsi 
%  wl
%  wh

% Data validation
% 
if real(z)<=0
 error('This is not a zero in the RHP')
end
%
if max(p)==0
  np=0;
else
       if min(real(p))<=0
            error('There is at least one pole not in the RHP')
       end
  np=length(p);
end
%
%Blaschke product for NMP poles
%
Bp=1;
 if np~=0
    for i=1:np
      Bp=abs(Bp*(p(i)-z)/(p(i)+z));
    end
 end
%
Ozl=atan((wl-imag(z))/real(z))+atan((wl+imag(z))/real(z));
Ozh=atan((wh-imag(z))/real(z))+atan((wh+imag(z))/real(z));
L1z=-log(epsi)*Ozl-log(1+epsi)*(pi-Ozh);
L2z=(-pi*log(Bp)+L1z)/(Ozh-Ozl);
sm=exp(L2z);







 

