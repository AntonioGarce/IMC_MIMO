function [a,b,c,d] = z2del(a,b,c,d,delta)
%Z2DEL   Conversion from discrete time z domain to discrete time delta domain.
%        [A,B,C,D] = Z2DEL(A,B,C,D,delta)
%        Given state space model in z domain returns state space maodel in 
%        delta domain.
%
%        [NUM,DEN] = Z2DEL(NUM,DEN,delta)
%        Given transfer function in z returns the equivalent tranfer function
%        in delta.

%        I.Kaspura
%        11/89


nargs = nargin;
if nargs == 3
    delta = c;
    [a,b,c,d] = tf2ss(a,b);
    nargs = 5;
end
error(nargchk(5,5,nargs));
error(abcdchk(a,b,c,d));

a = (a-eye(size(a)))/delta;
b = b/delta;

if nargout == 2
    [a,b] = ss2tf(a,b,c,d,1);
end


