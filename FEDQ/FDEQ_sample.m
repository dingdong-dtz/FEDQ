%This script is used to accurately calculate the balance distribution starting from rough EFIT balance
clear variables
close all
%
eq=init_from_grid(floor(64*1.5*2),floor(64*1.5*2*3));
paras=eq.paras;
paras.imax=15;


err_0=zeros(1,paras.imax);% Magnetic axis position error value
err_X=zeros(1,paras.imax);% Position error value of point X

delta_npsi=zeros(1,paras.imax);% Normalized magnetic flux variation
delta_psi=zeros(1,paras.imax);% Changes in absolute magnetic flux
delta_grid=zeros(1,paras.imax);% Changes in Grid Position
delta_0=zeros(1,paras.imax);% Changes in magnetic axis position
delta_X=zeros(1,paras.imax);% Change in the position of point X
%All spatial positions are dimensionless with large radii
eq.err=0.000000001;

icheck=0;
iit=1;
j=1;
nd=4;% (nd+1) is the number of lattice points used in the difference scheme for constructing the system of equations, and an even number is generally used to ensure central symmetry
md=9;% MD is the number of points used to calculate the metric difference format, usually using odd numbers to ensure center symmetry
order=4;% MD is the accuracy order of the differential format for calculating metric differences, usually using odd numbers to ensure center symmetry
X_n_fit=5;

eq.M = Metric_SOL_steady(eq.R,eq.Z,eq.rho,eq.theta,eq.geometry,md+4,order-2,X_n_fit,0);%Solving the metric
plot_matrix(eq.R,eq.Z,eq.psi,'psi',  eq.geometry,eq.boundary);
while (j<=paras.imax)&&iit==1
    %------------------ prvious status
    npsik=eq.npsi;  psik=eq.psi;
    Rk=eq.R; Zk=eq.Z; Rk0=eq.R0;Zk0=eq.Z0;Rkx=eq.Rx;Zkx=eq.Zx;
    %------------------- solve GS equation
    eq.Js = evalJs(eq.gnpsi,eq.Pprime,eq.FFprime,eq.gnpsi,eq.Pprime_pri,eq.FFprime_pri,eq.R,eq.npsi,eq.geometry);
    eq = Gsolver_SOL_speed(eq,nd);
    eq = upRZ_SOL_old(eq,icheck);
    eq.M = Metric_SOL_steady(eq.R,eq.Z,eq.rho,eq.theta,eq.geometry,md,order,X_n_fit,icheck);%Solving the metric
    %%
    %Find the relative change between two steps
    %Calculate the position changes of the magnetic axis and X point, as well as the maximum change in grid points
    dR=(eq.R-Rk)./eq.R0;    dZ=(eq.Z-Zk)./eq.R0;    dr=sqrt(dR.^2+dZ.^2);%Changes in grid points
    delta_grid(j)=max_geometry(dr,eq.geometry);
    dR0=eq.R0-Rk0;          dZ0=eq.Z0-Zk0;          dr0=sqrt(dR0.^2+dZ0.^2);%Changes in magnetic axis
    delta_0(j)=dr0./eq.R0;
    dRx=eq.Rx-Rkx;          dZx=eq.Zx-Zkx;          drx=sqrt(dRx.^2+dZx.^2);%Changes in point X
    delta_X(j)=drx./eq.R0;

    dpsi=(eq.psi-psik)/eq.psib;
    delta_psi(j)= max_geometry(abs(dpsi),eq.geometry);
    dnpsi=eq.npsi-npsik;
    delta_npsi(j)=max_geometry(abs(dnpsi),eq.geometry);

    disp('\Delta\psi_n')
    disp([num2str(delta_npsi(1:j),'%10.1e')])
    %%
    %------------------ check convergence
    if (delta_npsi(j)<eq.err)&&(delta_grid(j)<eq.err)
        iit=0;  %  finished
    else
        iit=1;  %  back to iteration
    end
    j=j+1;
end
%%
plot_contour(eq.R,eq.Z,eq.Rx,eq.Zx,eq.geometry,eq.boundary)
number=1:(j-1);
%%
%This section of code is used to generate the graphics required for the article
figure
pre=number;
semilogy(number(pre),abs(delta_npsi(pre)),'ko-',number(pre),delta_X(pre),'r+-',number(pre),delta_0(pre),'bx-','LineWidth',2);
legend('max($\delta\psi_n$)','$\delta\vec{r_x}$','$\delta\vec{r_0}$','Interpreter','latex')
xlabel('iteration step')
ylabel('Changes during iteration')
fff=gca;
kuan=fff.Position(3);
gao=fff.Position(4);
kuanfig=6.5;
gaofig=kuanfig/kuan*gao;
set(gcf,'unit','centimeters','Position',[12,2,kuanfig+1.5+0.5,gaofig+1.5+0.5]);  
set(gca,'unit','centimeters','Position',[1.5 1.5 kuanfig gaofig]); 
%exportgraphics(gcf,'convergence.png','Resolution',1000)

%%
function info=datainfo
%*************************************************
%  data information in the output structure eq
%*************************************************
info=struct('paras','a structure contains inputs parameter' ...
    ,'rho','normalized rho'...
    ,'theta','poloidal coordinate theta'...
    ,'R','solution of R(rho,theta)'...
    ,'Z','solution of Z(rho,theta)'...
    ,'R0','major radius of magnetic axis'...
    ,'Z0','Z position of magnetic axis'...
    ,'Rx','major radius of x point'...
    ,'Zx','Z position of x point'...
    ,'M','a structure contains metrics for the coordinates'...
    ,'Js','Js(rho,theta)=-[R^2 * mu0 * Pprime +  FFprime]'...
    ,'L','Coefficients for the GS operators'...
    ,'psi','Solution of GS equation psi_p'...
    ,'npsi','normalized psi_p with 0 at center and 1 at boundary'...
    ,'psib','psi_lcs - psi_axis'...
    ,'psi_lcs','psi at boudary'...
    ,'psi_axis','poloidal flux at magnetic axis'...
    ,'psi_x','poloidal flux at x point'...
    ,'psi_pri','poloidal flux at private region boudary'...
    ,'A','A matrix for A*x=B'...
    ,'B','B vector for A*x=B'...
    ,'nrho_axis','error in nrho at magnetic axis'...
    ,'boundary','storing boundary conditions'...
    ,'geometry','Data storage format'...
    ,'eq.r_range','Special processing range for singular points'...
    ,'gnpsi','Collection point of current pressure profile'...
    ,'eq.FFprime','Current and current gradient'...
    ,'eq.Pprime','pressure gradient'...
    ,'eq.gpsi_pri','Collection point of current pressure profile'...
    ,'eq.FFprime_pri','Current and current gradient'...
    ,'eq.Pprime_pri','pressure gradient'...
    );
end