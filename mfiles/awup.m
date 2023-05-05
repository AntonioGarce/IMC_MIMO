function [ci,nci,dci]=awup(nc,dc)
%
% FUNCTION [ci,nci,dci]=awup(nc,dc)
%
% Function to decompose a BIPROPER controller to implement
% an anti wind-up architecture.
%
% nc is the controller numerator polynomial
% dc is the controller denominator polynomial
% ci is the high frequency gain, i.e. ci=nc(1)/dc(1)
% nci and dci are polynomials which satisfy
%
%     nc      ci*dci
%    ---- =----------------
%     dc    dci  +   ci*nci
%



% First, eliminate any leading zeros
%
[nt,num]=p_elcero(nc);[dt,den]=p_elcero(dc);
%
% Check that the controller is biproper
%
if length(num)~=length(den)
  error ('Your controller is not biproper')
end
%
% High frequency gain
%
ci=num(1)/den(1);
%
% Computation of nci and dci
%
ncy=ci*den-num; dcy=ci*num;
ncyy=ncy/dcy(1); dci=dcy/dcy(1);
[nt,nci]=p_elcero(ncyy);

