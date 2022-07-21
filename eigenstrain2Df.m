function [V,D,P]=eigenstrain2Df(epsxx,epsyy,epsxy)
mat(1,1)=epsxx;
mat(1,2)=epsxy;
mat(2,1)=epsxy;
mat(2,2)=epsyy;
[Vtemp,Dtemp]=eig(mat);
V=[Vtemp(1,1),Vtemp(2,1) Vtemp(1,2),Vtemp(2,2)];
% Vtemp(1,1)*Vtemp(1,2)+Vtemp(2,1)*Vtemp(2,2)
% pause
D=[Dtemp(1,1) Dtemp(2,2)];
P=mean(D);
