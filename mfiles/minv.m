function [Ai,Bi,Ci,Di]=minv(A,B,C,D)
%
% FUNCTION [Ai,Bi,Ci,Di]=minv(A,B,C,D)
%
% This function allows to obtain the inverse of a 
% BIPROPER  MIMO system in state space form.
%
% See MATLAB commands linmod2 and dlinmod

% Is the given system biproper?
%
  if det(D)==0
    error('The given system is not BIPROPER')
  end
%
% State space form of the inverse
%
   Din=inv(D);
   Ai=A-B*Din*C;
   Bi=B*Din;
   Ci=-Din*C;
   Di=Din;


