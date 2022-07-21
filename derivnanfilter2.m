function [Ix,Iy] = derivnanfilter2(Ixe,Iye,pasxy,d)
%%% Attention il faut absolument que d soit normalis? correctement pour
%%% exemple voir gradientN_h2

iptsetpref('UseIPPL',false)% because of a matlab bug

if length(pasxy)==2
    pasx=pasxy(1);
    pasy=pasxy(2);
else
    pasx=pasxy;
    pasy=pasxy;
end
    
tab= floor((length(d)-1)/2);
Ix1=padarray(Ixe,[1 1],NaN);
Iy1=padarray(Iye,[1 1],NaN);


Ix=imfilter(Ix1/pasx,d,NaN);
Iy=imfilter(Iy1/pasy,d',NaN);


Ip1=padarray(Ix1,[tab tab],NaN);
Ixp=padarray(Ix,[tab tab],NaN);

Ip2=padarray(Iy1,[tab tab],NaN);
Iyp=padarray(Iy,[tab tab],NaN);


A=imfilter (Iy1,[1 -1]',NaN);
B=imfilter (Ix1,[1 -1],NaN);


dys=flipdim(Iy1,1);
dys2=flipdim(Ix1,2);
A2=imfilter (dys,[1 -1]');
B2=imfilter (dys2,[1 -1]);

switch length(d)
        case 2
        ind=isnan(A)&~isnan(Iy1);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind2(1:end-1,:)=ind(2:end,:);
        Iyp(ind)=(Ip2(ind2)-Ip2(ind))/pasxy;
        clear ind
        

        
        
        ind=isnan(B)&~isnan(Ix1);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind2(:,1:end-1)=ind(:,2:end);
        Ixp(ind)=(Ip1(ind2)-Ip1(ind))/pasxy;
        clear ind
        
        
        Iy=Iyp(tab+2:end-tab-1,tab+2:end-tab-1);
        Ix=Ixp(tab+2:end-tab-1,tab+2:end-tab-1);
  
        
    case 3
        
        ind=isnan(A)&~isnan(Iy1);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind2(1:end-1,:)=ind(2:end,:);
        Iyp(ind)=(Ip2(ind2)-Ip2(ind))/pasy;
        clear ind
        
        ind=isnan(A2)&~isnan(dys);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind2(1:end-1,:)=ind(2:end,:);
        ind=flipdim(ind,1);
        ind2=flipdim(ind2,1);
        Iyp(ind)=(Ip2(ind)-Ip2(ind2))/pasy;
        clear ind
        
        
        ind=isnan(B)&~isnan(Ix1);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind2(:,1:end-1)=ind(:,2:end);
        Ixp(ind)=(Ip1(ind2)-Ip1(ind))/pasx;
        clear ind
        
        ind=isnan(B2)&~isnan(dys2);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind2(:,1:end-1)=ind(:,2:end);
        ind=flipdim(ind,2);
        ind2=flipdim(ind2,2);
        Ixp(ind)=(Ip1(ind)-Ip1(ind2))/pasx;
        
        Iy=Iyp(tab+2:end-tab-1,tab+2:end-tab-1);
        Ix=Ixp(tab+2:end-tab-1,tab+2:end-tab-1);

        
        
    case 5
        C=imfilter(A,[1 -1]');
        D=imfilter(B,[1 -1]);

        C2=imfilter(A2,[1 -1]');
        D2=imfilter(B2,[1 -1]);
        %% premiere ligne de nan
        ind=isnan(A)&~isnan(Iy1);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind3=false(size(ind));
        ind2(1:end-1,:)=ind(2:end,:);
        ind3(1:end-2,:)=ind(3:end,:);
        Iyp(ind)=(d(2)*2*(Ip2(ind2)-Ip2(ind))+d(1)*2*(Ip2(ind3)-Ip2(ind)))/pasy;
        clear ind
        
        ind=isnan(A2)&~isnan(dys);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind3=false(size(ind));
        ind2(1:end-1,:)=ind(2:end,:);
        ind3(1:end-2,:)=ind(3:end,:);
        ind=flipdim(ind,1);
        ind2=flipdim(ind2,1);
        ind3=flipdim(ind3,1);
        Iyp(ind)=(d(2)*2*(Ip2(ind)-Ip2(ind2))+d(1)*2*(Ip2(ind)-Ip2(ind3)))/pasy;
        clear ind
        
        ind=isnan(B)&~isnan(Ix1);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind3=false(size(ind));
        ind2(:,1:end-1)=ind(:,2:end);
        ind3(:,1:end-2)=ind(:,3:end);
        Ixp(ind)=(d(2)*2*(Ip1(ind2)-Ip1(ind))+d(1)*2*(Ip1(ind3)-Ip1(ind)))/pasx;
        clear ind
        
        ind=isnan(B2)&~isnan(dys2);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind3=false(size(ind));
        ind2(:,1:end-1)=ind(:,2:end);
        ind3(:,1:end-2)=ind(:,3:end);
        ind=flipdim(ind,2);
        ind2=flipdim(ind2,2);
        ind3=flipdim(ind3,2);
        Ixp(ind)=(d(2)*2*(Ip1(ind)-Ip1(ind2))+d(1)*2*(Ip1(ind)-Ip1(ind3)))/pasx;
        
        %% deuxieme ligne de nan
        r=2*(d(1)+d(2))/(2*d(2)+d(1));
         
        ind=isnan(C)&~isnan(A);
        ind=padarray(ind,[tab tab]);
        
        ind2=false(size(ind));
        ind3=ind2;
        ind4=ind2;
        ind2(1:end-1,:)=ind(2:end,:);
        ind3(1:end-2,:)=ind(3:end,:);
        ind4(2:end,:)=ind(1:end-1,:);
        Iyp(ind)=(d(2)*r*(Ip2(ind2)-Ip2(ind4))+d(1)*r*(Ip2(ind3)-Ip2(ind)))/pasy;
        clear ind
        ind=isnan(C2)&~isnan(A2);
        ind=padarray(ind,[tab,tab]);
        ind2=false(size(ind));
        ind3=ind2;
        ind4=ind2;
        ind2(1:end-1,:)=ind(2:end,:);
        ind3(1:end-2,:)=ind(3:end,:);
        ind4(2:end,:)=ind(1:end-1,:);
        ind=flipdim(ind,1);
        ind2=flipdim(ind2,1);
        ind3=flipdim(ind3,1);
        ind4=flipdim(ind4,1);
        Iyp(ind)=(d(2)*r*(Ip2(ind4)-Ip2(ind2))+d(1)*r*(Ip2(ind)-Ip2(ind3)))/pasy;
        clear ind
        
        ind=isnan(D)&~isnan(B);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind3=ind2;
        ind4=ind2;
        ind2(:,1:end-1)=ind(:,2:end);
        ind3(:,1:end-2)=ind(:,3:end);
        ind4(:,2:end)=ind(:,1:end-1);
        Ixp(ind)=(d(2)*r*(Ip1(ind2)-Ip1(ind4))+d(1)*r*(Ip1(ind3)-Ip1(ind)))/pasx;       
        clear ind
        
        ind=isnan(D2)&~isnan(B2);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind3=ind2;
        ind4=ind2;
        ind2(:,1:end-1)=ind(:,2:end);
        ind3(:,1:end-2)=ind(:,3:end);
        ind4(:,2:end)=ind(:,1:end-1);
        ind=flipdim(ind,2);
        ind2=flipdim(ind2,2);
        ind3=flipdim(ind3,2);
        ind4=flipdim(ind4,2);
        Ixp(ind)=(d(2)*r*(Ip1(ind4)-Ip1(ind2))+d(1)*r*(Ip1(ind)-Ip1(ind3)))/pasx;
        
           Iy=Iyp(tab+2:end-tab-1,tab+2:end-tab-1);
         Ix=Ixp(tab+2:end-tab-1,tab+2:end-tab-1);
        
        
    case 7
        C=imfilter(A,[1 -1]');
        D=imfilter(B,[1 -1]);
        C2=imfilter(A2,[1 -1]');
        D2=imfilter(B2,[1 -1]);
        E=imfilter(C,[1 -1]');
        F=imfilter(D,[1 -1]);   
        E2=imfilter(C2,[1 -1]');
        F2=imfilter(D2,[1 -1]);
        
        %% premiere ligne de nan
        
        ind=isnan(A)&~isnan(Iy1);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind3=ind2;
        ind4=ind2;
        ind2(1:end-1,:)=ind(2:end,:);
        ind3(1:end-2,:)=ind(3:end,:);
        ind4(1:end-3,:)=ind(4:end,:);
        Iyp(ind)=(d(2)*2*(Ip2(ind3)-Ip2(ind))+d(1)*2*(Ip2(ind4)-Ip2(ind))+d(3)*2*(Ip2(ind2)-Ip2(ind)))/pasy;
        clear ind
        
        
        ind=isnan(A2)&~isnan(dys);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind3=ind2;
        ind4=ind2;
        ind2(1:end-1,:)=ind(2:end,:);
        ind3(1:end-2,:)=ind(3:end,:);
        ind4(1:end-3,:)=ind(4:end,:);
        ind=flipdim(ind,1);
        ind2=flipdim(ind2,1);
        ind3=flipdim(ind3,1);
        ind4=flipdim(ind4,1);
        Iyp(ind)=(d(2)*2*(Ip2(ind)-Ip2(ind3))+d(1)*2*(Ip2(ind)-Ip2(ind4))+d(3)*2*(Ip2(ind)-Ip2(ind2)))/pasy;
        clear ind
        
        ind=isnan(B)&~isnan(Ix1);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind3=ind2;
        ind4=ind2;
        ind2(:,1:end-1)=ind(:,2:end);
        ind3(:,1:end-2)=ind(:,3:end);
        ind4(:,1:end-3)=ind(:,4:end);
        Ixp(ind)=(d(2)*2*(Ip1(ind3)-Ip1(ind))+d(1)*2*(Ip1(ind4)-Ip1(ind))+d(3)*2*(Ip1(ind2)-Ip1(ind)))/pasx;
        clear ind
        
        ind=isnan(B2)&~isnan(dys2);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind3=ind2;
        ind4=ind2;
        ind2(:,1:end-1)=ind(:,2:end);
        ind3(:,1:end-2)=ind(:,3:end);
        ind4(:,1:end-3)=ind(:,4:end);
        ind=flipdim(ind,2);
        ind2=flipdim(ind2,2);
        ind3=flipdim(ind3,2);
        ind4=flipdim(ind4,2);
        Ixp(ind)=(d(2)*2*(Ip1(ind)-Ip1(ind3))+d(1)*2*(Ip1(ind)-Ip1(ind4))+d(3)*2*(Ip1(ind)-Ip1(ind2)))/pasx;
        
        %% deuxieme ligne de nan
        r=2*(d(1)+d(2)+d(3))/(2*d(3)+d(2));
        
        ind=isnan(C)&~isnan(A);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind3=ind2;
        ind4=ind2;
        ind5=ind2;
        ind2(1:end-1,:)=ind(2:end,:);
        ind3(1:end-2,:)=ind(3:end,:);
        ind4(1:end-3,:)=ind(4:end,:); 
        ind5(2:end,:)=ind(1:end-1,:); 
        Iyp(ind)=(d(3)*r*(Ip2(ind2)-Ip2(ind5))+d(1)*r*(Ip2(ind4)-Ip2(ind))+d(2)*r*(Ip2(ind3)-Ip2(ind)))/pasy;
        clear ind
        
        ind=isnan(C2)&~isnan(A2);
        ind=padarray(ind,[tab,tab]);
        ind2=false(size(ind));
        ind3=ind2;
        ind4=ind2;
        ind5=ind2;
        ind2(1:end-1,:)=ind(2:end,:);
        ind3(1:end-2,:)=ind(3:end,:);
        ind4(1:end-3,:)=ind(4:end,:);
        ind5(2:end,:)=ind(1:end-1,:);
        ind=flipdim(ind,1);
        ind2=flipdim(ind2,1);
        ind3=flipdim(ind3,1);
        ind4=flipdim(ind4,1);
        ind5=flipdim(ind5,1);
        Iyp(ind)=(d(3)*r*(Ip2(ind5)-Ip2(ind2))+d(1)*r*(Ip2(ind)-Ip2(ind4))+d(2)*r*(Ip2(ind)-Ip2(ind3)))/pasy;
        clear ind
        
        
        ind=isnan(D)&~isnan(B);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind3=ind2;
        ind4=ind2;
        ind5=ind2;
        ind2(:,1:end-1)=ind(:,2:end);
        ind3(:,1:end-2)=ind(:,3:end);
        ind4(:,1:end-3)=ind(:,4:end);
        ind5(:,2:end)=ind(:,1:end-1);
        Ixp(ind)=(d(3)*r*(Ip1(ind2)-Ip1(ind5))+d(1)*r*(Ip1(ind4)-Ip1(ind))+d(2)*r*(Ip1(ind3)-Ip1(ind)))/pasx;    
        clear ind
        
        ind=isnan(D2)&~isnan(B2);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind3=ind2;
        ind4=ind2;
        ind5=ind2;
        ind2(:,1:end-1)=ind(:,2:end);
        ind3(:,1:end-2)=ind(:,3:end);
        ind4(:,1:end-3)=ind(:,4:end);
        ind5(:,2:end)=ind(:,1:end-1);
        ind=flipdim(ind,2);
        ind2=flipdim(ind2,2);
        ind3=flipdim(ind3,2);
        ind4=flipdim(ind4,2);
        ind5=flipdim(ind5,2);
        Ixp(ind)=(d(3)*r*(Ip1(ind5)-Ip1(ind2))+d(1)*r*(Ip1(ind)-Ip1(ind4))+d(2)*r*(Ip1(ind)-Ip1(ind3)))/pasx;
               
        %% troisi?me ligne de nan
        r=2*(d(1)+d(2)+d(3))/(2*d(3)+2*d(2)+d(1));
        
        ind=isnan(E)&~isnan(C);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind3=ind2;
        ind4=ind2;
        ind5=ind2;
        ind6=ind2;
        ind2(1:end-1,:)=ind(2:end,:);
        ind3(1:end-2,:)=ind(3:end,:);
        ind4(1:end-3,:)=ind(4:end,:); 
        ind5(2:end,:)=ind(1:end-1,:); 
        ind6(3:end,:)=ind(1:end-2,:); 
        Iyp(ind)=(d(3)*r*(Ip2(ind2)-Ip2(ind5))+d(2)*r*(Ip2(ind3)-Ip2(ind6))+d(1)*r*(Ip2(ind4)-Ip2(ind)))/pasy;
        clear ind
        
        ind=isnan(E2)&~isnan(C2);
        ind=padarray(ind,[tab,tab]);
        ind2=false(size(ind));
        ind3=ind2;
        ind4=ind2;
        ind5=ind2;
        ind6=ind2;
        ind2(1:end-1,:)=ind(2:end,:);
        ind3(1:end-2,:)=ind(3:end,:);
        ind4(1:end-3,:)=ind(4:end,:);
        ind5(2:end,:)=ind(1:end-1,:);
        ind6(3:end,:)=ind(1:end-2,:);
        ind=flipdim(ind,1);
        ind2=flipdim(ind2,1);
        ind3=flipdim(ind3,1);
        ind4=flipdim(ind4,1);
        ind5=flipdim(ind5,1);
        ind6=flipdim(ind6,1);
        Iyp(ind)=(d(3)*r*(Ip2(ind5)-Ip2(ind2))+d(2)*r*(Ip2(ind6)-Ip2(ind3))+d(1)*r*(Ip2(ind)-Ip2(ind4)))/pasy;
        clear ind
        
        
        ind=isnan(F)&~isnan(D);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind3=ind2;
        ind4=ind2;
        ind5=ind2;
        ind6=ind2;
        ind2(:,1:end-1)=ind(:,2:end);
        ind3(:,1:end-2)=ind(:,3:end);
        ind4(:,1:end-3)=ind(:,4:end);
        ind5(:,2:end)=ind(:,1:end-1);
        ind6(:,3:end)=ind(:,1:end-2);
        Ixp(ind)=(d(3)*r*(Ip1(ind2)-Ip1(ind5))+d(2)*r*(Ip1(ind3)-Ip1(ind6))+d(1)*r*(Ip1(ind4)-Ip1(ind)))/pasx;
          clear ind
        
        ind=isnan(F2)&~isnan(D2);
        ind=padarray(ind,[tab tab]);
        ind2=false(size(ind));
        ind3=ind2;
        ind4=ind2;
        ind5=ind2;
        ind6=ind2;
        ind2(:,1:end-1)=ind(:,2:end);
        ind3(:,1:end-2)=ind(:,3:end);
        ind4(:,1:end-3)=ind(:,4:end);
        ind5(:,2:end)=ind(:,1:end-1);
        ind6(:,3:end)=ind(:,1:end-2);
        ind=flipdim(ind,2);
        ind2=flipdim(ind2,2);
        ind3=flipdim(ind3,2);
        ind4=flipdim(ind4,2);
        ind5=flipdim(ind5,2);
        ind6=flipdim(ind6,2);
        Ixp(ind)=(d(3)*r*(Ip1(ind5)-Ip1(ind2))+d(2)*r*(Ip1(ind6)-Ip1(ind3))+d(1)*r*(Ip1(ind)-Ip1(ind4)))/pasx;
                                  
        Iy=Iyp(tab+2:end-tab-1,tab+2:end-tab-1);
        Ix=Ixp(tab+2:end-tab-1,tab+2:end-tab-1);
end