function [sigxx,sigyy,sigxy]=straintostress2D(epsxx,epsyy,epsxy)

global E nu

factor1=E/(1-nu^2);
factor2=E/(1+nu);
sigxx=factor1*(epsxx+nu.*epsyy);
sigyy=factor1*(epsyy+nu.*epsxx);
sigxy=factor2.*epsxy;
