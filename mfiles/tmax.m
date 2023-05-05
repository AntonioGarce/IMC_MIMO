  function tm=tmax(p,z,epsi,wl,wh,tau)
%
% FUNCTION tm=tmax(p,z,epsi,wl,wh,tau)
%
% Function to compute a lower bound for the
% peak of the  nominal complementary sensitity To.
%
%  This is done in the following setup:
%
%  1)The plant model has a time delay tau,
%    a number of  (NMP) zeros in the open RHP,
%    and the effect of one  
%    pole in the open RHP is examined. 
%
%  2)Design specifications require that
%
%    |So|<epsi  for omega<wl
%    |To|<epsi  for omega >wh>wl
%
%  The user must provide the following data:
%
%  p : RHP  pole 
%  z : vector with NMP zeros (z=0 if no NMP zeros exists)
%  epsi
%  wl
%  wh
%
% Checking data validity
%
if real(p)<=0
 error('This is not a pole in the RHP')
end
%
if z==0
  nz=0;
else
     if min(real(z))<=0
          error('There is at least one zero not in the RHP')
     end
  nz=length(z);
end
%
%
%Blaschke product for NMP zeros
%
Bz=1;
 if nz~=0
    for i=1:nz
      Bz=abs(Bz*(z(i)-p)/(z(i)+p));
    end
 end
%
Opl=atan((wl-imag(p))/real(p))+atan((wl+imag(p))/real(p));
Oph=atan((wh-imag(p))/real(p))+atan((wh+imag(p))/real(p));
L1p=-log(epsi)*(pi-Oph)-log(1+epsi)*(Opl)+tau*real(p);
%
L2p=(-pi*log(Bz)+L1p)/(Oph-Opl);
tm=exp(L2p);







 

