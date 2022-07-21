clc
close all;
clear ;

%% Define some constants

eps1=0.00000001;

nu=0.49; %Poisson modulus
pix=1;% pixel size in µm


E=1000; %Young modulus


param.E=E;
param.nu=nu;
param.pix=pix;
param.eps1=eps1;



alpha = 0.33; % Size of arrow head relative to the length of the vector
beta = 0.33;  % Width of the base of the arrow head relative to the length 


load('cell_rect_2gaussians_truncated_ad_unplan.txt')
x=cell_rect_2gaussians_truncated_ad_unplan(:,1);
y=cell_rect_2gaussians_truncated_ad_unplan(:,2);
fact=1e6;% to convert loaded displacement data which are in meter in micrometer
dx=cell_rect_2gaussians_truncated_ad_unplan(:,3)*fact;
dy=cell_rect_2gaussians_truncated_ad_unplan(:,4)*fact;

Ny=100;
Nx=100;
crp=25;
indc=x>-crp&x<crp&y>-crp&y<crp;

%% fmult to have the result for a force of 1kPa the simulation was run with 1Pa but everything is linear!
fmult=1000;

figure
quiver(x,y,dx,dy,'g')

axis equal

N=2^9;
param.N=N;

Lx=pix*crp;
Ly=pix*crp;




Fx=scatteredInterpolant(x(indc),y(indc),dx(indc));
Fy=scatteredInterpolant(x(indc),y(indc),dy(indc));

pas=2*crp/N;
[X2,Y2]=meshgrid(-crp:pas:(N-1)/N*crp,-crp:pas:(N-1)/N*crp);


%%
UX=Fx(X2,Y2);
UY=Fy(X2,Y2);

dimx=size(UX,2);
dimy=size(UY,1);

Fx.Method = 'nearest';
Fy.Method= 'nearest';
UX(isnan(UX)) = Fx(X2(isnan(UX)), Y2(isnan(UX)));
UY(isnan(UY)) = Fy(X2(isnan(UY)), Y2(isnan(UY)));
UX=UX*pix;
UY=UY*pix;

[px,sx]=perdecomp(UX);
[py,sy]=perdecomp(UY);




Y=fft2(px,size(X2,1),size(X2,2));
Z=fft2(py,size(X2,1),size(X2,2));


Y=fftshift(Y);
Z=fftshift(Z);

kx=(-(size(X2,2))*pi+eps1:2*pi:(size(X2,2)-1)*pi)/(Lx);
ky=(-(size(X2,1))*pi+eps1:2*pi:(size(X2,1)-1)*pi)/(Ly);



[KX,KY]=meshgrid(kx,ky);

A=2*(1+nu)/E;


K=sqrt(KX.^2+KY.^2);

C=1./(A*(1-nu).*K);
K1=C.*((1-nu).*K.^2+nu*KX.^2);
K2=C.*nu.*KX.*KY;
K3=C.*((1-nu)*K.^2+nu*KY.^2);


    VX=K1.*Y+K2.*Z;
    VY=K2.*Y+K3.*Z;

VX=ifftshift(VX);
VY=ifftshift(VY);


ResultX(:,:)=real(ifft2(VX));
ResultY(:,:)=real(ifft2(VY));


F=sqrt(ResultX.^2+ResultY.^2);

figure('Name','Forces_gaussian')
    figToolbarFix
surf(X2,Y2,F*fmult,'edgecolor','none')
shading interp
colorbar;
view(0,90)
    axis tight
        axis equal
set(gca,'Xlim',[-12,12],'Ylim',[-12 12],'Clim',[0 100])



%%
figure('Name','quivergaussad')
pas=5;
hold on
figToolbarFix
[~,autoscalef]=quiverh(X2(1:pas:end,1:pas:end),Y2(1:pas:end,1:pas:end),ResultX(1:pas:end,1:pas:end)*fmult,ResultY(1:pas:end,1:pas:end)*fmult,'g');
axis equal;
axis tight;

   set(gca,'Ydir','reverse')

set(gca,'Xlim',[-12 12],'Ylim',[-12 12])
x0=-10;
y0=-11;
alpha=0.33;
beta=0.33;
maxim2=500;
u=maxim2*autoscalef;
v=0;
uu = [x0;x0+u];
vv = [y0;y0+v];
plot(uu,vv,'k','LineWidth',2);
  hu = [x0+u-alpha*(u+beta*(v+eps));x0+u; ...
          x0+u-alpha*(u-beta*(v+eps))];
  hv = [y0+v-alpha*(v-beta*(u+eps));y0+v; ...
        y0+v-alpha*(v+beta*(u+eps))];
plot(hu,hv,'k','LineWidth',2);
text(x0+u+2,y0,[num2str(maxim2),'  Pa'],'BackgroundColor',[1 1 1]);

%%

save(['force_fft_cell_rect_ad_2gaussians_',num2str(N),'.mat'], 'X2','Y2','ResultX','ResultY','param')



