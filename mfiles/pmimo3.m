%This program obtain the Q matrix for the solved problem 13.2
%
%The user should provide D_11 (s) (relative degree 1)
%and  D_22 (s) (relative degree 2)
%
nd11=input('Enter numerator for D_{11}(s):   '); 
dd11=input('Enter denominator for D_{11}(s):   ');
%
nd22=input('Enter numerator for D_{22}(s):   '); 
dd22=input('Enter denominator for D_{22}(s):   ');
%
nd21=conv([1 -5],nd11);dd21=conv([1 5],dd11);
%
nq11=conv([1 3],nd21);dq11=dd21;
%
nq12=conv([1 3],nd22);dq12=dd22;
%
nq21=10*conv([1 4 3],nd11);dq21=conv([1 5],dd11);
%
nq22=-conv([1 4 3],nd22);dq22=dd22;

