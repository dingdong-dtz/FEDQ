function [eq,error] = FDEQ_theory(eq)
%   This function is used to construct the solution validation of analytical equilibrium
%   This function omits many result displays
%%
%The following program is used to determine the absolute accurate position of the magnetic axis and X point using analytical methods
%Firstly, find the position of point X by seeking a local minimum
[R0_real,Z0_real,psi0_theory,Rx_real,Zx_real,~] = find_real_axis_x(eq.R0,eq.Z0,eq.Rx,eq.Zx);
%%
paras=eq.paras;
paras.imax=15;%Set the expected number of iterations
abs_err_psi_mean=zeros(1,paras.imax);
abs_err_psi_var=zeros(1,paras.imax);
abs_err_psi_Max=zeros(1,paras.imax);
err_0=zeros(1,paras.imax);
err_X=zeros(1,paras.imax);

delta_psi=zeros(1,paras.imax);
delta_grid=zeros(1,paras.imax);
delta_0=zeros(1,paras.imax);
delta_X=zeros(1,paras.imax);
%All spatial positions are dimensionless with large radii
eq.err=1e-15;

icheck=0;
iit=1;
j=1;
nd=4;%（nd+1）It is the number of lattice points used in the difference scheme for constructing the system of equations, 
% and an even number is generally used to ensure central symmetry
md=9;%mdIt is the number of points used to calculate the differential format of the metric, usually using odd numbers to ensure center symmetry
order=4;%mdIt is the accuracy order of the differential format when calculating the metric, usually using odd numbers to ensure center symmetry
X_n_fit=5;
while (j<=paras.imax)&&iit==1
    %------------------ prvious status
    npsik=eq.npsi;  psik=eq.psi;
    Rk=eq.R; Zk=eq.Z; Rk0=eq.R0;Zk0=eq.Z0;Rkx=eq.Rx;Zkx=eq.Zx;%These are used to determine the difference between two steps to determine if the loop can be stopped
    %------------------- solve GS equation
    %The normalized magnetic flux current pressure function should be selected, 
    % otherwise it faces extrapolation problems and does not conform to the physical meaning of the profile
    eq.Js = evalJs(eq.gnpsi,eq.Pprime,eq.FFprime,eq.gnpsi,eq.Pprime_pri,eq.FFprime_pri,eq.R,eq.npsi,eq.geometry);
    eq.M = Metric_SOL_steady(eq.R,eq.Z,eq.rho,eq.theta,eq.geometry,md,order,X_n_fit,icheck);%Solving the metric
    eq = Gsolver_SOL_speed(eq,nd);
    %%
    %------------------- update [R,Z] grids
    %eq.psi=theory_sample(eq.R,eq.Z);
    eq = upRZ_SOL_old(eq,icheck);
    %%
    %Verify the absolute error after solving each step of the program
    psi_theory = theory_sample(eq.R,eq.Z);
    error=abs(eq.psi-psi_theory)/eq.psib;%Normalize the absolute error
    abs_err_Mmax = max_geometry(error,eq.geometry);%This item represents the maximum absolute value error
    abs_err_psi_Max(j)=abs_err_Mmax;

    abs_err_Mmean = mean_geometry(error,eq.geometry,eq.M.J,eq.R,eq.Z,eq.rho,eq.theta);%This item represents the average absolute value error
    abs_err_psi_mean(j)=abs_err_Mmean;

    abs_err_var = sqrt(mean_geometry(error.^2,eq.geometry,eq.M.J,eq.R,eq.Z,eq.rho,eq.theta));%This item represents the root mean square error
    abs_err_psi_var(j)=abs_err_var;
    %Calculate the position error between the magnetic axis and point X, and use large radius normalization
    err_R0=eq.R0-R0_real;          err_Z0=eq.Z0-Z0_real;          err_r0=sqrt(err_R0.^2+err_Z0.^2);%Error situation of magnetic axis
    err_0(j)=err_r0./eq.R0;
    err_Rx=eq.Rx-Rx_real;          err_Zx=eq.Zx-Zx_real;          err_rx=sqrt(err_Rx.^2+err_Zx.^2);%Error situation of point X
    err_X(j)=err_rx./eq.R0;
    %%
    %Find the relative change between two steps
    dpsi=(eq.psi-psik)/eq.psib;
    delta_psi(j)= max_geometry(abs(dpsi),eq.geometry);

    %Calculate the position changes of the magnetic axis and X point, as well as the maximum change in grid points
    dR=(eq.R-Rk)./eq.R0;    dZ=(eq.Z-Zk)./eq.R0;    dr=sqrt(dR.^2+dZ.^2);%Changes in grid points
    delta_grid(j)=max_geometry(dr,eq.geometry);
    dR0=eq.R0-Rk0;          dZ0=eq.Z0-Zk0;          dr0=sqrt(dR0.^2+dZ0.^2);%Changes in magnetic axis
    delta_0(j)=dr0./eq.R0;
    dRx=eq.Rx-Rkx;          dZx=eq.Zx-Zkx;          drx=sqrt(dRx.^2+dZx.^2);%Changes in point X
    delta_X(j)=drx./eq.R0;
    %%
    j=j+1;
end
%Save the workspace for solving
save([num2str(eq.geometry.nr) '_' num2str(eq.geometry.nt)]);
error=zeros(6,1);
error(1)=floor(sqrt(eq.geometry.nr*(eq.geometry.t_max-eq.geometry.t_min+1)));%Actual grid point density
error(2)=abs(abs_err_psi_Max(end));
error(3)=abs(abs_err_psi_mean(end));
error(4)=abs(abs_err_psi_var(end));
error(5)=err_0(end);
error(6)=err_X(end);
error(7)=eq.r_range;

end