%Function to eliminate the leading zeros from a given vector vfat.
%The systen returns the lean vector vec and its dimension  ntrue
%
function [ntrue,vec]=p_elcero(vfat);
[x,y]=size(vfat);
ntrue=y;
i=1;
while ((abs(vfat(i))<1e-5)&(i<(y+1))), 
   ntrue=ntrue-1;
   i=i+1;
end;

for k=1:ntrue,
  vec(k)=vfat(y-ntrue+k);
end;
