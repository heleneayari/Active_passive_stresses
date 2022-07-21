%% Code of BISM/ KISM (Vincent NIER, 15/07/16)
%From a movie of traction force data, input: TFMfield.mat, this program
%calculates the internal stress tensor, output: stress.mat.
%System:
%-R, length on y direction between ymin and ymax
%-C, length on x direction between xmin and xmax
%-kM, number of time step
%Different algorithms: BISM, KISM (Kalman condition) and different
%hyperparameters calculation methods: MAP, L-curve, ... (meth_lambda condition,
%meth_omega) are possible. Stress boundary conditions can be imposed (BC
%condition).
%If you use BISM/ KISM please cite the following reference in your work:
%Nier, V. Inference of internal stress in a cell monolayer, Biophysical
%Journal, 110(7), 1625-1635, 2016


clear;
close all;
clc

r=groot;
sz=r.ScreenSize;
pos=[10 10 sz(4),sz(4)];

pix=1;%µm/pixel
Nforce=512;
% nrect=1;
%%






% system size for data analysis

kM=1;% number of timestep


%% Parameters
BC =0;
% 1: free stress boundary conditions \sigma_{ik}n_k=0 imposed in the prior
% 0: no boundary conditions imposed in the prior
meth_lambda = 3;
stepLambda=3;
% method to determine the regularization parameter Lambda
% 1: MAP
% 2: L-curve
% 3: fixed value
if meth_lambda == 3
    fixed_lambda =0.4;  % lambda value when fixed
end

noise_value= 0.;
% noise amplitude of the traction force data
% (same unit as traction force data)
% set to zero when unknown
fplot = 0;
% Graphics
% 0: no figures
% 1: plot all figures
% 2: plot only stress related figures
% 3: plot stress tensor figure
mult=1;
% multiplicative coefficient for a better graphical representation of the stress tensor
stepM=10;
%number of step for the iteration of the hyperparameters, MAP method
Kalman=0; %0: BISM, 1: KISM
meth_omega=2; %method to determine the regularization parameter \Omega if
%Kalman=1, 1: determination from traction force variation, 2: set to a fixed
%value




load(['force_fft_cell_rect_ad_2gaussians_',num2str(Nforce),'.mat']);
% define spatial grid
coeff=pix; %conversion  coefficient pixel to microns
Lcell=10;
xmin=-coeff*Lcell;
xmax=coeff*Lcell;
ymin=-coeff*Lcell;
ymax=coeff*Lcell;
% spatial resolution (assume same resolution in x and y
% deltax=deltay=delta)
C = 51; % number of columns

delta =(xmax-xmin)/(C);


gx=xmin:delta:xmax;
gy=ymin:delta:ymax;
gx2=xmin-delta:delta:xmax+delta;
gy2=xmin-delta:delta:xmax+delta;
R=length(gy)-1;
N = R*C; %system size for traction force field
N2=R*C+2*R+2*C; %size when adding boundaries conditions
Ninf=4*N+2*(C+R); %system size for stress tensor

x=(gx(1:end-1)+gx(2:end))/2;
y=(gy(1:end-1)+gy(2:end))/2;

%for fixed_lambda=[ 1e-4 1e-3 2e-3 5e-3 1e-2 1e-1 1]
for fixed_lambda=[3e-3]
    fixed_lambda
    
    %% Matrices A and B
    
    
    %% Building 2N*Ninf matrix A such that A*\sigma=T (A discretized form of the
    %    divergence operator)
    %Building matrix Ax (derivative wrt x)
    Ax1=diag(ones(C-1, 1), 1)-diag(ones(C, 1), 0); %bulk (order 2)
    Ax1=[Ax1,[zeros(C-1, 1); 1]];
    
    Ax=sparse(Ax1);
    for i=1:R-1
        Ax=blkdiag(Ax,Ax1); %iteration L times for all the rows
    end
    clear Ax1;
    %Building matrix Ay (derivative wrt y)
    Ay1=diag(ones((R-1)*C,1), C)-diag(ones(R*C,1 ), 0); %bulk (order 2)
    Ay=sparse([Ay1, [zeros((R-1)*C,C); eye(C)]]);
    clear Ay1;
    
    
    Bbcx=zeros(2*C, C*(R+1));
    Bbcy=zeros(2*R, (C+1)*R);
    Bbcx(1:C,1:C)=speye(C);
    Bbcx(C+1:2*C, 1+R*C:(R+1)*C)=speye(C);
    for i=1:R
        Bbcy(i, 1+(i-1)*(C+1))=1;
        Bbcy(i+R, C+1+(i-1)*(C+1))=1;
    end
    Bbc=sparse([Bbcy, sparse(2*R,(C+1)*R+2*C*(R+1)); sparse(2*C,(C+1)*R),Bbcx, sparse(2*C,(C+1)*R+C*(R+1)); sparse(2*C,(C+1)*R+C*(R+1)), Bbcx,sparse(2*C,(C+1)*R); sparse(2*R,(C+1)*R+2*C*(R+1)), Bbcy]); %xx, yy, yx and xy
    %Building A from Ax and Ay
    A1=[Ax, sparse(N,N+C) Ay, sparse(N,N+R); sparse(N,N+R) Ay, sparse(N,N+C),Ax]/delta;
    A=[A1;Bbc];
    clear Ax Ay;
    %% Building Ninf*Ninf B covariance matrix of the prior, S_0^{-1}=B/s_0^2
    %% Stress norm regularization
    B=speye(Ninf);
    %% Add prior enforcing the equality of shear components
    d1=diag(ones((R-1)*C, 1), C)+diag(ones(R*C,1), 0); %bulk (order 2)
    bd1=[d1,[zeros((R-1)*C, C); eye(C)]];
    clear d1;
    d2=-diag(ones(C-1, 1), 1)-diag(ones(C, 1),0); %bulk (order 2)
    d2=[d2,[zeros(C-1, 1); -1]];
    bd2=d2;
    for i=1:R-1
        bd2=blkdiag(bd2, d2); %iteration L times for all the rows
    end
    clear d2;
    Bdiff=sparse([zeros(N, 2*N+C+R), bd1, bd2]);
    clear bd1 bd2;
    alpha_xy=10^3; %hyperparameter on the shear equality (has to be large)
    B=B+alpha_xy^2*sparse(Bdiff'*Bdiff);
    clear Bdiff;
    %% Add prior enforcing the boundary conditions
    if BC==1
        Bbcx=zeros(2*C, C*(R+1));
        Bbcy=zeros(2*R, (C+1)*R);
        Bbcx(1:C,1:C)=speye(C);
        Bbcx(C+1:2*C, 1+R*C:(R+1)*C)=speye(C);
        for i=1:R
            Bbcy(i, 1+(i-1)*(C+1))=1;
            Bbcy(i+R, C+1+(i-1)*(C+1))=1;
        end
        Bbc=sparse([Bbcy, sparse(2*R,(C+1)*R+2*C*(R+1)); sparse(2*C,(C+1)*R),Bbcx, sparse(2*C,(C+1)*R+C*(R+1)); sparse(2*C,(C+1)*R+C*(R+1)), Bbcx,sparse(2*C,(C+1)*R); sparse(2*R,(C+1)*R+2*C*(R+1)), Bbcy]); %xx, yy, yx and xy
        %components
        clear Bbcx Bbcy;
        alpha_BC=10^3; %hyperparameter on the boundary condition (has to be large)
        B=B+alpha_BC^2*sparse(Bbc'*Bbc);
        clear Bbc
    end
    
    %% Traction force
    for k0=1:kM %iteration on the time
        load(['force_fft_cell_rect_ad_2gaussians_',num2str(Nforce),'.mat']);
        
        % average the force on the BISM grid
        [Nb,edgesx,edgesy,BINx,BINy]=histcounts2(X2*pix,Y2*pix,gx2,gy2);
        for uu=1:size(Nb,1)
            for vv=1:size(Nb,2)
                ind=BINx==uu&BINy==vv;
                mTx(vv,uu)=mean(ResultX(ind));
                mTy(vv,uu)=mean(ResultY(ind));
                
            end
        end

        vTx=reshape(mTx(2:end-1,2:end-1)', N, 1); %Tx in a N vector form
        vTy=reshape(mTy(2:end-1,2:end-1)', N, 1); %Ty in a N vector form
        T=[vTx; vTy]; %vector T (2N components)
        
        %% We will now add the boundaries conditions:
        Tb=[-(mTx(2:end-1,1)+mTx(2:end-1,2))/2;(mTx(2:end-1,C+2)+mTx(2:end-1,C+1))/2;-(mTy(1,2:end-1)'+mTy(2,2:end-1)')/2;(mTy(R+1,2:end-1)'+mTy(R+2,2:end-1)')/2;-(mTy(2:end-1,1)+mTy(2:end-1,2))/2;(mTy(2:end-1,C+2)+mTy(2:end-1,C+1))/2;-(mTx(1,2:end-1)'+mTx(2,2:end-1)')/2;(mTx(R+2,2:end-1)'+mTx(R+1,2:end-1)')/2];
        
        
        T=[T;-Tb];
        %% Hyperparameters determination methods
        %% Lambda
        if noise_value==0
            s=1; %fixed to constant value
        else
            s=noise_value; %settled to the experimental value
        end
        if meth_lambda==1
            %MAP estimation with Jeffreys hyperprior
            stepM=5; %number of iterative steps
            % initiate iteration
            step=1;
            beta(step)=10^(0); %beta=s^2
            alpha(step)=delta^2*beta(1)/10^(-4); %alpha=s0^2
            while step<stepM
                sigmaMAP=(beta(step)*B/alpha(step)+(A'*A))\(A'*T);
                sigmaMAP_norm(step)=sigmaMAP'*B*sigmaMAP;
                res_norm(step)=(T-A*sigmaMAP)'*(T-A*sigmaMAP);
                step=step+1;
                alpha(step)=sigmaMAP_norm(step-1)/(4*N+2*(C+L)+2);
                beta(step)=res_norm(step-1)/(2*N+2);
            end
            figure
            semilogy(0:stepM-1,(beta./alpha).^0.5,'LineWidth',2)
            set(gca, 'FontSize', 18, 'fontName','Times');
            set(gcf,'Color','w')
            xlabel('Step','Fontsize',18)
            ylabel('\lambda (\mum^{-1})', 'Fontsize',18)
            s0=alpha(stepM)^(1/2);
            s=beta(stepM)^(1/2);
        elseif meth_lambda==2
            %L-curve
            k=1;
            for l=-6:0.1:1
                lambda_Lcurve(k)=10^l;
                Lsigma_post=(lambda_Lcurve(k)^2*B+(A'*A))\(A'*T);
                sigma_norm(k)=Lsigma_post'*B*Lsigma_post;
                res_norm(k)=(T-A*Lsigma_post)'*(T-A*Lsigma_post);
                k=k+1;
            end
            dL=(log(sigma_norm(3:end))-log(sigma_norm(1:end-2)))./(2*log(res_norm(2:end-1))); %L curve slope
            ddL=(log(sigma_norm(3:end))-2*log(sigma_norm(2:end-1))+log(sigma_norm(1:end-2)))./(log(res_norm(2:end-1)).^2); %l curve curvature
            figure %plot of the L-curve
            loglog(res_norm/2,sigma_norm/2,'LineWidth',2)
            set(gca, 'FontSize', 18, 'fontName','Times');
            set(gcf,'Color','w')
            xlabel('$\mathrm{Residual\,norm\,\Vert\vec{T}-\textbf{A}\vec{\sigma}\Vert\,(kPa)}$','interpreter','latex','Fontsize',18)
            ylabel('Prior norm(kPa.\mum)', 'Fontsize',18)
            figure %plot of the L-curve slope evolution
            semilogx(res_norm(2:end-1)/2,dL,'LineWidth',2)
            set(gca, 'FontSize', 18, 'fontName','Times');
            set(gcf,'Color','w')
            xlabel('$\mathrm{Residual\,norm\,\Vert\vec{T}-\textbf{A}\vec{\sigma}\Vert\,(kPa)}$','interpreter','latex','Fontsize',18)
            ylabel('$\mathrm{\frac{d\,log(\vec{\sigma}^T\textbf{B}\vec{\sigma})^{1/2}}{d\,log\Vert\vec{T}-\textbf{A}\vec{\sigma}\Vert} (\mu m)}$','interpreter','latex','Fontsize',18)
            figure %plot of the L-curve curvature evolution
            semilogx(res_norm(2:end-1)/2,ddL,'LineWidth',2)
            set(gca, 'FontSize', 18, 'fontName','Times');
            set(gcf,'Color','w')
            xlabel('$\mathrm{Residual\,norm\,\Vert\vec{T}-\textbf{A}\vec{\sigma}\Vert\,(kPa)}$','interpreter','latex','Fontsize',18)
            ylabel('$\mathrm{\frac{d^2\,log(\vec{\sigma}^T\textbf{B}\vec{\sigma})^{1/2}}{d\,log^2\Vert\vec{T}-\textbf{A}\vec{\sigma}\Vert} (\mu m.kPa^{-1})}$','interpreter','latex','Fontsize',18)
            [max_curv,index_max]=max(ddL);
            lambda=lambda_Lcurve(index_max+1); %lambda defined at the maximum of the L-curve curvature
            s0=s/lambda;
        elseif meth_lambda==3
            %Fix Lambda value
            Lambda=fixed_lambda;
            s0=s*delta/Lambda^0.5;
        end
        %% Omega
        if Kalman==1
            if meth_omega==1
                var_tractionx=mean(mean((mTx- TFMfield.( sprintf('frame%d',k0+1).Tx).^2)));
                var_tractiony=mean(mean((mTy- TFMfield.( sprintf('frame%d',k0+1).Ty).^2)));
                var_traction=(var_tractionx+var_traction_y)/2;
                Omega=var_traction/(4*s^2)*Lambda;
            elseif meth_omega==2
                Omega=10^2*Lambda;
            end
            gamma=s0*Omega^0.5;
        end
        %% Save the hyperparameters
        stress_BISM.s{k0}=s;
        stress_BISM.s0{k0}=s0;
        if Kalman==1
            stress_BISM.gamma{k0}=gamma;
        end
        %% Covariance definition
        S_0inv=sparse(B)/s0^2; %Covariance matrix Sigma_0^{-1}=B/s0^2
        Sinv=speye(2*N2)/s^2; %Covariance noise matrix Sigma{-1}
        S=s^2*speye(2*N2);
        if (Kalman==1) && (k0>1)
            Cinv=speye(Ninf)+gamma^2*S_0inv;
            sigma_post_evo=Cinv\sigma_post_obs{k0-1};
            Spost_evo=sparse((Cinv\Spost_obs+gamma^2*speye(Ninf))*(Cinv\speye(Ninf)));
            iSpost_evo=sparse((Spost_obs*(Cinv\speye(Ninf))+gamma^2*speye(Ninf))\Cinv);
            sigma_post_obs{k0}=sigma_post_evo+(s^2*iSpost_evo+A'*A)\(A'*(T-A*sigma_post_evo));
            Spost_obs=sparse((speye(Ninf)-((s^2*iSpost_evo+A'*A)\(A'*A)))*Spost_evo);
            clear Spost_evo iSpost_evo sigma_post_evo;
        else
            Spost_obs=sparse((S_0inv+A'*Sinv*A)\speye(Ninf));
            sigma_post_obs{k0}=Spost_obs*(A'*Sinv*T);
        end
        if noise_value>0
            Spost{k0}=Spost_obs;
        end
    end
    clear B Sinv S S_0inv;
    %% Data interpolation and stress representation
    for k0=1:kM
        %Extract xx, yy, xy and yx components
        if noise_value>0
            Spost_xx=reshape(diag(Spost{k0}(1:N+R, 1:N+R))', C+1, R)'; %inferred
            variance on stress xx component
            Spost_yy=reshape(diag(Spost{k0}(N+R+1:2*N+C+R, N+R+1:2*N+C+R))', C, R+1)';
            %inferred variance on stress yy component
            Spost_xy=reshape(diag(Spost{k0}(2*N+C+R+1:3*N+2*C+R,2*N+C+R+1:3*N+2*C+R))', C, R+1)'; %inferred variance on stress xy component
            Spost_yx=reshape(diag(Spost{k0}(3*N+2*C+R+1:4*N+2*(C+R),3*N+2*C+R+1:4*N+2*(C+R)))', C+1, R)'; %inferred variance on stress xy component
        end
        vsigma_post_xx=sigma_post_obs{k0}(1:N+R);
        vsigma_post_yy=sigma_post_obs{k0}(N+R+1:2*N+C+R);
        vsigma_post_xy=sigma_post_obs{k0}(2*N+C+R+1:3*N+2*C+R);
        vsigma_post_yx=sigma_post_obs{k0}(3*N+2*C+R+1:4*N+2*(C+R));
        %Interpolation on the grid of the force data
        for i=1:R
            for j=1:C
                vsigma_xx((i-1)*C+j, 1)=(vsigma_post_xx((i-1)*(C+1)+j)+vsigma_post_xx((i-1)*(C+1)+j+1))/2;
                vsigma_yx((i-1)*C+j, 1)=(vsigma_post_yx((i-1)*(C+1)+j)+vsigma_post_yx((i-1)*(C+1)+j+1))/2;
                vsigma_yy((i-1)*C+j, 1)=(vsigma_post_yy((i-1)*C+j)+vsigma_post_yy((i-1)*C+j+C))/2;
                vsigma_xy((i-1)*C+j, 1)=(vsigma_post_xy((i-1)*C+j)+vsigma_post_xy((i-1)*C+j+C))/2;
                if noise_value>0
                    Ssigma_xx(i,j)=(Spost_xx(i, j)+Spost_xx(i, j+1)+2*Spost((i-1)*(C+1)+j, (i-1)*(C+1)+j+1))/4;
                    Ssigma_yx(i,j)=(Spost_yx(i, j)+Spost_yx(i,j+1)+2*Spost(3*N+2*C+R+(i-1)*(C+1)+j, 3*N+2*C+R+(i-1)*(C+1)+j+1))/4;
                    Ssigma_yy(i,j)=(Spost_yy(i, j)+Spost_yy(i+1, j)+2*Spost(N+R+(i-1)*C+j, N+R+i*C+j))/4;
                    Ssigma_xy(i,j)=(Spost_xy(i, j)+Spost_xy(i+1,j)+2*Spost(2*N+C+R+(i-1)*C+j, 2*N+C+R+i*C+j))/4;
                end
            end
        end
        vsigma_post_xx=vsigma_xx;
        vsigma_post_yy=vsigma_yy;
        vsigma_post_xy=(vsigma_yx+vsigma_xy)/2;
        if noise_value>0
            Spost_xx=Ssigma_xx;
            Spost_yy=Ssigma_yy;
            Spost_xy=(Ssigma_xy+Ssigma_yx)/4;
            clear Ssigma_xx Ssigma_yy Ssigma_xy Ssigma_yx Spost;
            spostxx=reshape(Spost_xx.^0.5, N, 1);
            spostyy=reshape(Spost_yy.^0.5, N, 1);
            spostxy=reshape(Spost_xy.^0.5, N, 1);
        end
        %Reshape the stress vector in a matrix form
        sigma_post_xx=reshape(vsigma_post_xx, C, R)'; %inferred stress xx
        sigma_post_yy=reshape(vsigma_post_yy, C, R)'; %inferred stress yy
        sigma_post_xy=reshape(vsigma_post_xy, C, R)'; %inferred stress xy
        %% Inferred stress tensor figure
        if fplot==1
            for i=1:R
                for j=1:C
                    [V,D]=eig([sigma_post_xx(i,j) sigma_post_xy(i,j); sigma_post_xy(i,j) sigma_post_yy(i,j)]); %eigenvalue and eigenvector calculus
                    ab1=[x(j)-mult*D(1,1)*V(1,1),x(j)+mult*D(1,1)*V(1,1)];
                    or1=[y(i)-mult*D(1,1)*V(2,1),y(i)+mult*D(1,1)*V(2,1)];
                    if D(1,1)>=0
                        figure(k0)
                        hold on
                        plot(ab1,or1,'b','LineWidth',2);
                    else
                        figure(k0)
                        hold on
                        plot(ab1,or1,'r','LineWidth',2);
                    end
                    ab2=[x(j)-mult*D(2,2)*V(2,1),x(j)+mult*D(2,2)*V(2,1)];
                    or2=[y(i)+mult*D(2,2)*V(1,1),y(i)-mult*D(2,2)*V(1,1)];
                    if D(2,2)>=0
                        figure(k0)
                        hold on
                        plot(ab2,or2,'b','LineWidth',2);
                    else
                        figure(k0)
                        hold on
                        plot(ab2,or2,'r','LineWidth',2);
                    end
                end
            end
            set(gca, 'FontSize', 18, 'fontName','Times');
            set(gcf,'Color','w')
            xlabel('x (\mum)', 'Fontsize', 18)
            ylabel('y (\mum)', 'Fontsize', 18)
            title('\sigma^{inf}', 'Fontsize', 18)
            axis([xmin-0.01*(xmax-xmin) xmax+0.01*(xmax-xmin) ymin-0.01*(ymax-ymin)
                ymax+0.01*(ymax-ymin)])
        end
        %% Quantitative evaluation of inference
        Tinf=A*sigma_post_obs{k0}; %Calculate Tinf =A*\tau_inf=A*tpost
        clear sigma_post_obs{k0};
        Tinfx=Tinf(1:N); Tinfy=Tinf(N+1:2*N);

        
        R2x=1-((vTx'-Tinfx')*(vTx-Tinfx))/((vTx'-mean(vTx))*(vTx-mean(vTx)));
        R2y=1-((vTy'-Tinfy')*(vTy-Tinfy))/((vTy'-mean(vTy))*(vTy-mean(vTy)));
        R2=(R2x+R2y)/2; %determination coefficient between the data and the traction
        %obtained from the BISM stress (should be close to 1)
        clear Tinf;
        %Stress mean values from traction force
        if BC==1
            sxx_mean=-sum(sum(reshape(vTx, C, R)'.*(ones(R, 1)*(x-(xmax-xmin)/2))))/N;
            syy_mean=-sum(sum(reshape(vTy, C, R)'.*((y-(ymax-ymin)/2)'*ones(1,C))))/N;
            sxy_mean=-sum(sum(reshape(vTy, C, R)'.*(ones(R, 1)*(x-(xmax-xmin)/2))+reshape(vTx, C, R)'.*((y-(ymax-ymin)/2)'*ones(1,C))))/(2*N);
            stress.mean_from_t{k0}=[sxx_mean,syy_mean,sxy_mean];
        end
    end
    
    stress.syy=sigma_post_yy;
    stress.sxx=sigma_post_xx;
    stress.sxy=sigma_post_xy;
    stress.x=x;
    stress.y=y;
    stress.R2_T=R2;
    save(['Stress_BISM_ad_2gaussians_lambda_',num2str(Lambda), '_cart_BCp.mat'] ,'stress');
    
end



