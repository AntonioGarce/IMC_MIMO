disp('This programme generates all parameters to run newss.mdl ');
disp('The user must previously design biproper controllers for y and z');
P=input('Enter a row vector with  four observer poles:   ');
[nu,du,ny,dy]=lambor(P);
cyn=input('Enter the numerator of the Y biproper controler as a row vector:   ');
cyd=input('Enter the denominator of the Y biproper controler as a row vector:   ');
[cyi,nhy,dhy]=awup(cyn,cyd);
czn=input('Enter the numerator of the Z biproper controler as a row vector:   ');
czd=input('Enter the denominator of the Z biproper controler as a row vector:   ');
[czi,nhz,dhz]=awup(czn,czd);
zsp=input('Enter desired set point for state z:   ');
zh=input('Enter the upper limit for z:   ');
zl=input('Enter the lower limit for z:   ');
zsat=input('Enter saturation value for the magnitude of the state z:   ');
