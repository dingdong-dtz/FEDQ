clear variables
close all
%
%eq=initial_from_theory(floor(64*2^1.8),floor(64*4*2^1.8));
eq=initial_from_theory(floor(64),floor(64*4));
%%
%The following program is used to determine the absolute accurate position of the magnetic axis and X point using analytical methods, in order to calculate errors
%Firstly, the position of the magnetic axis is determined by finding the local minimum, and the position of point X is determined by finding the zero point of the derivative function
[R0_real, Z0_real,psi0_theory,Rx_real,Zx_real,psi_X_theory] = find_real_axis_x(eq.R0,eq.Z0,eq.Rx,eq.Zx);
%%
paras=eq.paras;
paras.imax=10;% Set the expected number of iterations
abs_err_psi_mean=zeros(1,paras.imax);% Mean error value
abs_err_psi_var=zeros(1,paras.imax);% Mean square root error value
abs_err_psi_Max=zeros(1,paras.imax);% Maximum error value
abs_nerr_psi_mean=zeros(1,paras.imax);% Mean error value
abs_nerr_psi_var=zeros(1,paras.imax);% Mean square root error value
abs_nerr_psi_Max=zeros(1,paras.imax);% Maximum error value
err_0=zeros(1,paras.imax);% Magnetic axis position error value
err_X=zeros(1,paras.imax);% Position error value of point X
delta_npsi=zeros(1,paras.imax);% Normalized magnetic flux variation
delta_psi=zeros(1,paras.imax);% Changes in absolute magnetic flux
delta_grid=zeros(1,paras.imax);% Changes in Grid Position
delta_0=zeros(1,paras.imax);% Changes in magnetic axis position
delta_X=zeros(1,paras.imax);% Change in the position of point X
%All spatial positions are dimensionless with large radii
eq.err=1e-15;
icheck=0;
iit=1;
j=1;
nd=4;% （nd+1）It is the number of lattice points used in the difference scheme for constructing the system of equations,
% and an even number is generally used to ensure central symmetry
md=9;% mdIt is the number of points used to calculate the differential format of the metric, usually using odd numbers to ensure center symmetry
order=4;% mdIt is the accuracy order of the differential format when calculating the metric, usually using odd numbers to ensure center symmetry
X_n_fit=5;
eq.M = Metric_SOL_steady(eq.R,eq.Z,eq.rho,eq.theta,eq.geometry,md,order,X_n_fit,icheck);% Solving the metric
while (j<=paras.imax)&&iit==1
%------------------ prvious status
npsik=eq.npsi;  psik=eq.psi;
Rk=eq.R; Zk=eq.Z; Rk0=eq.R0; Zk0=eq.Z0; Rkx=eq.Rx; Zkx=eq.Zx;% These are used to determine the difference between two steps to determine if the loop can be stopped
%------------------- solve GS equation
eq.Js = evalJs(eq.gnpsi,eq.Pprime,eq.FFprime,eq.gnpsi,eq.Pprime_pri,eq.FFprime_pri,eq.R,eq.npsi,eq.geometry);% Equation source term
eq = Gsolver_SOL_speed(eq,nd);% Solving difference equations
eq = upRZ_SOL_old(eq,icheck);% Update magnetic mesh
eq.M = Metric_SOL_steady(eq.R,eq.Z,eq.rho,eq.theta,eq.geometry,md,order,X_n_fit,icheck);% Solving the metric
%%
%Verify the absolute error after solving each step of the program
%Due to the use of fitting methods to update the magnetic axis and magnetic flux values at point X during the process of updating the grid, this may have a certain reduction effect on errors
psi_theory = theory_sample(eq.R,eq.Z);
npsi_theory=(psi_theory-psi0_theory)/(eq.psi_lcs-psi0_theory);
error=abs(eq.psi-psi_theory)/eq.psib;% Normalize the absolute error
nerror=abs(eq.npsi-npsi_theory);% Normalized magnetic flux error
abs_err_Mmax = max_geometry(error,eq.geometry);% This item represents the maximum absolute value error
abs_err_psi_Max(j)=abs_err_Mmax;
abs_err_Mmean = mean_geometry(error,eq.geometry,eq.M.J,eq.R,eq.Z,eq.rho,eq.theta);% This item represents the average absolute value error
abs_err_psi_mean(j)=abs_err_Mmean;
abs_err_var = sqrt(mean_geometry(error.^2,eq.geometry,eq.M.J,eq.R,eq.Z,eq.rho,eq.theta));% This item represents the root mean square error
abs_err_psi_var(j)=abs_err_var;
%
abs_nerr_Mmax = max_geometry(nerror,eq.geometry);% This item represents the maximum absolute value error
abs_nerr_psi_Max(j)=abs_nerr_Mmax;
abs_nerr_Mmean = mean_geometry(nerror,eq.geometry,eq.M.J,eq.R,eq.Z,eq.rho,eq.theta);% This item represents the average absolute value error
abs_nerr_psi_mean(j)=abs_nerr_Mmean;
abs_nerr_var = sqrt(mean_geometry(nerror.^2,eq.geometry,eq.M.J,eq.R,eq.Z,eq.rho,eq.theta));% This item represents the root mean square error
abs_nerr_psi_var(j)=abs_nerr_var;
%Calculate the position error between the magnetic axis and point X, and use large radius normalization
err_R0=eq.R0-R0_real;          err_Z0=eq.Z0-Z0_real;          err_r0=sqrt(err_R0.^2+err_Z0.^2);% Error situation of magnetic axis
err_0(j)=err_r0./eq.R0;
err_Rx=eq.Rx-Rx_real;          err_Zx=eq.Zx-Zx_real;          err_rx=sqrt(err_Rx.^2+err_Zx.^2);% Error situation of point X
err_X(j)=err_rx./eq.R0;
%%
%Find the relative change between two steps
%Calculate the position changes of the magnetic axis and X point, as well as the maximum change in grid points
dR=(eq.R-Rk)./eq.R0;    dZ=(eq.Z-Zk)./eq.R0;    dr=sqrt(dR.^2+dZ.^2);% Changes in grid points
delta_grid(j)=max_geometry(dr,eq.geometry);
dR0=eq.R0-Rk0;          dZ0=eq.Z0-Zk0;          dr0=sqrt(dR0.^2+dZ0.^2);% Changes in magnetic axis
delta_0(j)=dr0./eq.R0;
dRx=eq.Rx-Rkx;          dZx=eq.Zx-Zkx;          drx=sqrt(dRx.^2+dZx.^2);% Changes in point X
delta_X(j)=drx./eq.R0;
dpsi=(eq.psi-psik)/eq.psib;% Normalize and dimensionless the changes in \ psi itself, and compare% before the grid changes;
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
plot_contour(eq.R,eq.Z,eq.Rx,eq.Zx,eq.geometry,eq.boundary);% Draw a grid distribution map
%%
number=1:(j-1);
%Generate error convergence graph in the article
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
exportgraphics(gcf,'convergence.png','Resolution',1000)
%%
%Distribution of absolute error, illustrated in the article
psi_theory = theory_sample(eq.R,eq.Z);
error=abs(eq.psi-psi_theory)/eq.psib+realmin;
plot_matrix_delta(eq.R,eq.Z,log10(error),'delta',eq.geometry,eq.boundary,[-3:-0.25:-7]);
c=colorbar;
c.Label.String = 'Log(\psi_{err})'; c.Label.FontSize=10;
axis equal
boundary=eq.boundary;
kuan=max(max(boundary.Rlcs),max(boundary.Rpri))-min(min(boundary.Rlcs),min(boundary.Rpri));
gao=-min(min(boundary.Zlcs),min(boundary.Zpri))+max(max(boundary.Zlcs),max(boundary.Zpri));
xlim([min(min(boundary.Rlcs),min(boundary.Rpri))-kuan*0.05 max(max(boundary.Rlcs),max(boundary.Rpri))+kuan*0.05])
ylim([min(min(boundary.Zlcs),min(boundary.Zpri))-gao*0.05 max(max(boundary.Zlcs),max(boundary.Zpri))+gao*0.05])
kuanfig=4;
gaofig=kuanfig/kuan*gao;
set(gcf,'unit','centimeters','Position',[12,2,kuanfig+1+0.5+1.5,gaofig+1.2+0.5]);
set(gca,'unit','centimeters','Position',[1 1.2 kuanfig gaofig]);
exportgraphics(gcf,'error.png','Resolution',1000);
%%
%Draw an image to demonstrate the matching between magnetic surface and analytical solution
plot_contour_comparing(eq.R,eq.Z,eq.Rx,eq.Zx,eq.geometry,eq.boundary,eq.rho.^2*eq.psib+eq.psi_axis,psi_X_theory)
%%
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
,' psi','Solution of GS equation psi_p'...
,'npsi','normalized psi_p with 0 at center and 1 at boundary'...
,'psib','psi_lcs - psi_axis'...
,'psi_lcs','psi at boudary'...
,'psi_axis','poloidal flux at magnetic axis'...
,'psi_x','poloidal flux at x point'...
,'psi_pri','poloidal flux at private region boudary'...
,'A','A matrix for A*x=B'...
,' B','B vector for A*x=B'...
,' nrho_axis','error in nrho at magnetic axis'...
,'boundary','storing boundary conditions'...
,' geometry','Data storage format'...
,'eq.r_range','Special processing range for singular points'...
,'gnpsi','Collection point of current pressure profile'...
,'eq.FFprime','Current and current gradient'...
,'eq.Pprime','pressure gradient'...
,'eq.gpsi_pri','Collection point of current pressure profile'...
,' eq.FFprime_pri','Current and current gradient'...
,'eq.Pprime_pri','pressure gradient'...
);
end
