function [a,b,c,d] = del2z(a,b,c,d,delta)
%DEL2Z   Conversion from discrete time delta domain to discrete time z domain.
%        [A,B,C,D] = DEL2Z(A,B,C,D,delta)
%        Given state space model in delta domain and the sampling frequency, delta, returns the state 
%        space model in z domain.
%
%        [NUM,DEN] = DEL2Z(NUM,DEN,delta)
%        Given transfer function in delta domain and the sampling frequency, delta, returns the 
%        equivalent tranfer function in z. 

%        11/89
%        I.Kaspura

%        Further modifications by Rick Middleton, 8/91, to ensure it works
%        properly with first order state space systems.

nargs = nargin;
if nargs == 3
    delta = c;
    [a,b,c,d] = tf2ss(a,b);
    nargs = 5;
end
error(nargchk(5,5,nargs));
error(abcdchk(a,b,c,d));

n=size(a);
I=eye(n(1),n(1));

a = a*delta + I;
b = b*delta;

if nargout == 2
    [a,b] = ss2tf(a,b,c,d,1);
end

