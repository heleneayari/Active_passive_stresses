clear
close all;
clc;
global E nu;
E=5000;% value used for the rigidity of the cell in the comsol simulation
nu=0.49;


load('cell_rect_2gaussians_truncated_ad_unplan.txt')
x=cell_rect_2gaussians_truncated_ad_unplan(:,1);
y=cell_rect_2gaussians_truncated_ad_unplan(:,2);
% fact is used to convert displacement in micrometers they are in meters in
% the table
fact=1e6;
dxc=cell_rect_2gaussians_truncated_ad_unplan(:,3)*fact;
dyc=cell_rect_2gaussians_truncated_ad_unplan(:,4)*fact;

scrsz = get(0,'ScreenSize');



Mx=max(x);
mx=min(x);
My=max(y);
my=min(y);
k=1;
%for pasxy=[3 10 16 20 24 32 48 64]
for pasxy=[0.25]

[X,Y]=meshgrid([mx:pasxy:Mx],[my:pasxy:My]);
[Xc,Yc]=meshgrid([mx+pasxy/2:pasxy:Mx-pasxy/2],[mx+pasxy/2:pasxy:My-pasxy/2]);



    
    %% methode 1 : Interpolation
    %ind=1:length(res(1).sp);
    Fx=scatteredInterpolant(x,y,dxc);
    Fy=scatteredInterpolant(x,y,dyc);

     dxi=Fx(X,Y);
    dyi=Fy(X,Y);   


    
    
    
            [epsxx,epsxy]=gradientN_h2(dxi,'fb',pasxy);
            [epsyx,epsyy]=gradientN_h2(dyi,'fb',pasxy);


    epsxy=1/2*(epsxy+epsyx);
 
    
    [vec, Dt, P1] = arrayfun(@eigenstrain2Df,epsxx,epsyy,epsxy,'UniformOutput',0);
    AA1=cell2mat(vec(:));
    BB1=cell2mat(Dt(:));
    P=cell2mat(P1);

    
    [vecdev, Dtdev,~] = arrayfun(@eigenstrain2Df,epsxx-P,epsyy-P,epsxy,'UniformOutput',0);
    AA=cell2mat(vecdev(:));
    BB=cell2mat(Dtdev(:));
    eps(k).P=P(:);
    eps(k).X=Xc(:);
    eps(k).Y=Yc(:);
    
    eps(k).vecdevx1=AA(:,1);
    eps(k).vecdevy1=AA(:,2);
    eps(k).vecdevx2=AA(:,3);
    eps(k).vecdevy2=AA(:,4);
    
    eps(k).vecx1=AA1(:,1);
    eps(k).vecy1=AA1(:,2);
    eps(k).vecx2=AA1(:,3);
    eps(k).vecy2=AA1(:,4);
    eps(k).lambdamin=BB1(:,1);
    eps(k).lambdamax=BB1(:,2);
    eps(k).lambdadev=BB(:,2);
    
    [sigxx,sigyy,sigxy]=arrayfun(@straintostress2D,epsxx,epsyy,epsxy);
    [vec, Dt, P1] = arrayfun(@eigenstrain2Df,sigxx,sigyy,sigxy,'UniformOutput',0);
    AA1=cell2mat(vec(:));
    BB1=cell2mat(Dt(:));
    P=cell2mat(P1);
    [vecdev, Dtdev, P1dev] = arrayfun(@eigenstrain2Df,sigxx-P,sigyy-P,sigxy,'UniformOutput',0);
    AA=cell2mat(vecdev(:));
    BB=cell2mat(Dtdev(:));

    stress(k).sigxx=sigxx;
    stress(k).siggxxm=min(sigxx(:));
    stress(k).siggxxM=max(sigxx(:));
    stress(k).siggyym=min(sigyy(:));
    stress(k).siggyyM=max(sigyy(:));
    stress(k).siggxym=min(sigxy(:));
    stress(k).siggxyM=max(sigxy(:));
    stress(k).sigyy=sigyy;
    stress(k).sigxy=sigxy;
    stress(k).P=P(:);

    stress(k).X=X;
    stress(k).Y=Y;
    stress(k).vecdevx1=AA(:,1);
    stress(k).vecdevy1=AA(:,2);
    stress(k).vecdevx2=AA(:,3);
    stress(k).vecdevy2=AA(:,4);
    stress(k).vecx1=AA1(:,1);
    stress(k).vecy1=AA1(:,2);
    stress(k).vecx2=AA1(:,3);
    stress(k).vecy2=AA1(:,4);
    stress(k).lambdamin=BB1(:,1);
    stress(k).lambdamax=BB1(:,2);
    stress(k).lambdadev=BB(:,2);
    stress(k).nstress=sqrt(sigxx.^2+sigyy.^2+2*sigxy.^2);
    stress(k).nstressM=max((stress(k).nstress(:)));
    
    
    

save(['resultepsadrectad2gaussians_pas',num2str(pasxy),'_nu',num2str(nu),'.mat'],'eps')
save(['resultstressadrectad2gaussians_pas',num2str(pasxy),'_nu',num2str(nu),'.mat'],'stress')

end
