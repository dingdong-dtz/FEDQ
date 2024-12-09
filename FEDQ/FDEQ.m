function eq = FDEQ(eq)

paras=eq.paras;
paras.imax=15;

delta_npsi=zeros(1,paras.imax);
delta_psi=zeros(1,paras.imax);
delta_grid=zeros(1,paras.imax);
delta_0=zeros(1,paras.imax);
delta_X=zeros(1,paras.imax);
%All spatial positions are dimensionless with large radii
eq.err=0.00000000000000000000001;

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
    Rk=eq.R; Zk=eq.Z; Rk0=eq.R0;Zk0=eq.Z0;Rkx=eq.Rx;Zkx=eq.Zx;
    %------------------- solve GS equation
    eq.Js = evalJs(eq.gnpsi,eq.Pprime,eq.FFprime,eq.gnpsi,eq.Pprime_pri,eq.FFprime_pri,eq.R,eq.npsi,eq.geometry);
    eq.M = Metric_SOL_steady(eq.R,eq.Z,eq.rho,eq.theta,eq.geometry,md,order,X_n_fit,icheck);
    eq = Gsolver_SOL_speed(eq,nd);
    %%
    %------------------- update [R,Z] grids
    %eq.psi=theory_sample(eq.R,eq.Z);
    eq = upRZ_SOL_old(eq,icheck);
    %%
    %Find the relative change between two steps
    %Calculate the position changes of the magnetic axis and X point, as well as the maximum change in grid points
    dR=(eq.R-Rk)./eq.R0;    dZ=(eq.Z-Zk)./eq.R0;    dr=sqrt(dR.^2+dZ.^2);%网格点的变化情况
    delta_grid(j)=max_geometry(dr,eq.geometry);
    dR0=eq.R0-Rk0;          dZ0=eq.Z0-Zk0;          dr0=sqrt(dR0.^2+dZ0.^2);%磁轴的变化情况
    delta_0(j)=dr0./eq.R0;
    dRx=eq.Rx-Rkx;          dZx=eq.Zx-Zkx;          drx=sqrt(dRx.^2+dZx.^2);%X点的变化情况
    delta_X(j)=drx./eq.R0;

    dpsi=(eq.psi-psik)/eq.psib;
    delta_psi(j)= max_geometry(abs(dpsi),eq.geometry);
    dnpsi=eq.npsi-npsik;
    delta_npsi(j)=max_geometry(abs(dnpsi),eq.geometry);
    %%
    %------------------ check convergence
    if (delta_npsi(j)<eq.err)&&(delta_grid(j)<eq.err)
        %iit=0;  %  finished
    else
        iit=1;  %  back to iteration
    end
    j=j+1;
end
end