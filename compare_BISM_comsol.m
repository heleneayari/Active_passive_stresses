clc
close all
clear
%%

makePretty;

EE=1000;
pix=1;




name='Large_gauss_ad';
namestress='resultstressadrectad2gaussians_pas';
nameBISM='Stress_BISM_ad_2gaussians_l';
namecomsol='cell_rect_2gaussians_ad_truncated_stressplane.txt';
namecomsolcell='cell_alone_rect_2gaussians_truncated_stressplane.txt';
paramreg=[1];
fact=1e3;% Comsol simulation was run with a force of 1Pa, we make it here 1kPa as everything is linear
x=[0:0.1:10];
y=x;
[X,Y]=meshgrid(x,y);


    %% Load Comsol results cell+ substrate
    T=readtable(namecomsol);
    
    A=table2array(T);
    Fxx=scatteredInterpolant(A(:,1),A(:,2),A(:,3));
    Fyy=scatteredInterpolant(A(:,1),A(:,2),A(:,4));
    Fxy=scatteredInterpolant(A(:,1),A(:,2),A(:,5));
    sxx=Fxx(X,Y);
    syy=Fyy(X,Y);
    sxy=Fxy(X,Y);
    
    ff=figure(1000);
    ff.Name=[name,'FEM_ISM'];
    figToolbarFix
    hxxm=plot(x,sxx(1,:)*fact);
    

    hold on

    
    %% Load Comsol results cell+ alone
    clear T A Fxx Fyy Fxy
    T=readtable(namecomsolcell);
    [X,Y]=meshgrid([0:0.1:10],[0:0.1:10]);
    A=table2array(T);
    Fxx=scatteredInterpolant(A(:,1),A(:,2),A(:,3));
    Fyy=scatteredInterpolant(A(:,1),A(:,2),A(:,4));
    Fxy=scatteredInterpolant(A(:,1),A(:,2),A(:,5));
    sxx0=Fxx(X,Y);
    syy0=Fyy(X,Y);
    sxy0=Fxy(X,Y);
    


    %%

    ff2=figure(100);
    ff2.Name=[name,'FEM_BISM'];

    figToolbarFix
    hold on
    plot(x,(sxx(1,:)-sxx0(1,:))*fact,'k')    
    
    
    h=figure(500);
    h.Name='comsolwwogel';
    figToolbarFix
    hold on
    plot(x,sxx(1,:)*fact)
    plot(x,sxx0(1,:)*fact) 
    legend({'cell+gel','cell alone'})
    xlabel('x ($\mu m$)')
    ylabel('$(Sc)_{xx}$ (Pa)')
    r=get(gcf,'position');
    set(gcf,'position',[r(1)/1.4 r(2)/4 2*r(3) 2*r(4)])
    leg(1)={'FEM'};
    
    
    %% load BISM

    fich=dir([nameBISM,'*cart_BCp.mat']);

    nu=0.49;
    pasxy=0.25;
    
    %% Load stress ISM and calculate div
    
    si=load([namestress,num2str(pasxy),'_nu',num2str(nu),'.mat']);
    % let's interpol stress on BISM grid
    N=(length(si.stress.X)+1)/2;

    cp=2;
    %for uu=1:length(fich)
    for uu=paramreg
        
        uu
        load(fich(uu).name);

        deb=length(nameBISM)+7;
        lambda(uu)=str2double(fich(uu).name(deb:end-13));
 
   
        figure(100)
        plot(stress.x(26:end),stress.sxx(26,26:end)*fact)

        leg(cp)={num2str(lambda(uu))};
 
        cp=cp+1;
    end



    
  
    xx=linspace(0,10,100);
    sxxism=interp1(si.stress.X(N,:),si.stress.sigxx(N,:)*fact,xx);
    figure(1000)
    
    cut=2;
    hism=plot(xx(1:end-cut),sxxism(1:end-cut));
    leg2(1)={'FEM'};
    leg2(2)={'ISM'};
    ylabel('$(S_c)_{xx}$ (Pa)')
    xlabel('$x$ ($\mu m$)')
    leg2f=legend(leg2);
    leg2f.Position=[  0.2106    0.7304    0.2928    0.1520];
    set(gcf,'position',[r(1)/1.4 r(2)/4 2*r(3) 2*r(4)])


    

    figure(100)

    ylabel('$(S_{tot})_{xx}$ (Pa)')
    xlabel('$x$ ($\mu m$)')

    set(gca,'Ylim',[0 300])
    legend1=legend(leg);
    r=get(gcf,'position');
    set(gcf,'position',[r(1)/1.4 r(2)/4 2*r(3) 2*r(4)])
    set(legend1,...
        'Position',[0.708115854055173 0.159725573430435 0.223899719469835 0.374576282905319]);
    

    
    
    
    
 