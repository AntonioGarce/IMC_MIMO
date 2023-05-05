function [a,b,c,d] = c2del(a,b,c,d,delta)
%C2DEL     Conversion from continuous domain to discrete time delta domain.
%          [DA,DB,DC,DD] = C2DEL(A,B,C,D,delta)  converts the continuous
%          time state space system:
%		     .
%		     x = Ax + Bu
%            y = Cx + Du
%
%  	       to the discrete time state space system:
%
%		     px = DAx + DBu
%             y = DCx + DDu
%
%	       using delta as the sampling period.
%
%          [DNUM,DDEN] = C2DEL(NUM,DEN,delta)  converts the continous system
%          y(s)/u(s) = NUM(s)/DEN(s)  to the discrete time delta domain 
%          system  y(d)/u(d) = NUM(d)/DEN(d).


%	   I.Kaspura
%	   11/89


nargs = nargin;
if nargs == 3                      % If in transfer function form transform
    delta = c;                     % to state space form.
    [a,b,c,d] = tf2ss(a,b);
    nargs = 5;
end

error(nargchk(5,5,nargs));
error(abcdchk(a,b,c,d));

if delta ~= 0                      % If delta=0 then no change occurs. 
    [m,n] = size(a);
    [m,nb] = size(b);
    I = eye(n,n); O = 0*I;
    omega = [I O]*expm(([a I; O O])*delta)*[O;I]/delta;
    %  omega is integral from 0 to delta of exp(a t) dt  over delta
    a = omega*a;
    b = omega*b;
end
       
if nargout == 2                    % If output is to be in transfer function
    [a,b] = ss2tf(a,b,c,d,1);      % form transform back.
end

