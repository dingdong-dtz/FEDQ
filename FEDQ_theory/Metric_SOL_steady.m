function M=Metric_SOL_steady(R,Z,rho,theta,geometry,number,order,X_n_fit,icheck)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%M=Metric_SOL_steady(R,Z,rho,theta,geometry,number,order,X_n_fit,icheck)
%This function uses a robust difference scheme to calculate the metric of the magnetic surface coordinate system
%This function uses the results of analytical fitting to process the metric near singular points
%This function does not smooth the partial derivatives obtained
%The number of grid points used in the number difference format
%When calculating the metric, only the first-order derivative format was actually used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nt2=geometry.nt2;
nt=geometry.nt;
nr_inner=geometry.nr_inner;
nr_outer=geometry.nr_outer;
nr_down=geometry.nr_down;
nr=geometry.nr;
rho_x=rho(nr_inner+1);
t_min=geometry.t_min;
t_max=geometry.t_max;
tr_max=geometry.tr_max;
tl_min=geometry.tl_min;
[d1x,kdx,d1y,kdy]=diff_me_steady_change_speed(rho,theta,number,order,geometry);

err=1.0e-16;  %  change 0 to err
%   Calculation of partial derivatives in the core area
[dRdr,dRdt]=partial(R,kdx,kdy,d1x,d1y,geometry);
[dZdr,dZdt]=partial(Z,kdx,kdy,d1x,d1y,geometry);
%    plot_matrix_near_x(R,Z,(dZdt),'dZdt',geometry,[])

sep_range=2;
for_mean=nr_inner+1-sep_range:nr_inner+1+sep_range;
d0=diff_steady_all(rho_x,rho(for_mean),0,2);
for_change=[nt1+1-floor(number/2):nt1+1+floor(number/2) nt1+nt_inner+1-floor(number/2):nt1+nt_inner+1+floor(number/2)];
dRdt(nr_inner+1,for_change)=d0*dRdt(for_mean,for_change);
dZdt(nr_inner+1,for_change)=d0*dZdt(for_mean,for_change);
%    plot_matrix_near_x(R,Z,(dZdt),'dZdt',geometry,[])

dRdt(1,:)=0;dZdt(1,:)=0;

J =  R.*(dRdr.*dZdt-dRdt.*dZdr);    
J=Eliminating_Singularities(J,err);
Jni=1./J;

cJ  = J./R.^2;
grt0 = -(dRdt.*dRdr+dZdt.*dZdr);
gtt0 = (dRdr.^2+dZdr.^2);
grr0 = (dRdt.^2+dZdt.^2);
cJ1 = (R./J).^2;
grt = cJ1.*grt0;
gtt = cJ1.*gtt0;
grr = cJ1.*grr0;
%%

%Update the metric using analytical and fitting information
[gtt_th,dtdR,dtdZ] = gtt_theory(R,Z,theta,rho,geometry);
dd=floor(number/2);
gtt(1:nr,tr_max+dd+1:tl_min-dd-1)=gtt_th(1:nr,tr_max+dd+1:tl_min-dd-1);
[grr_th,drdR,drdZ,index_0,index_x,weight]=fit_axis_and_x(R,Z,rho,theta,geometry,X_n_fit);
grr(index_x)=weight(index_x).*(grr_th(index_x)-grr(index_x))+grr(index_x);
grr(index_0)=weight(index_0).*(grr_th(index_0)-grr(index_0))+grr(index_0);
%Invert Jacobian and update the results of magnetic axis and X
Jni_th=Jni;
Jni_th(index_x)=1./R(index_x).*(drdR(index_x).*dtdZ(index_x)-drdZ(index_x).*dtdR(index_x));
Jni_th(index_0)=1./R(index_0).*(drdR(index_0).*dtdZ(index_0)-drdZ(index_0).*dtdR(index_0));
Jni_th(1,nt1+1:nt1+nt_inner)=10^16;
Jni(index_x)=weight(index_x).*(Jni_th(index_x)-Jni(index_x))+Jni(index_x);
Jni(index_0)=weight(index_0).*(Jni_th(index_0)-Jni(index_0))+Jni(index_0);
JupR2=Jni.*R.^2;
%Update the magnetic axis and X-point position of the mixing degree gauge
grt_th=grt;
grt_th(index_0)=(drdR(index_0).*dtdR(index_0)+drdZ(index_0).*dtdZ(index_0));
grt_th(index_x)=(drdR(index_x).*dtdR(index_x)+drdZ(index_x).*dtdZ(index_x));
grt_th(1,nt1+1:nt1+nt_inner)=10^16;
grt(index_x)=weight(index_x).*(grt_th(index_x)-grt(index_x))+grt(index_x);
grt(index_0)=weight(index_0).*(grt_th(index_0)-grt(index_0))+grt(index_0);

%%
crr = grr0./J;
crr(index_x)=grr(index_x)./Jni(index_x)./R(index_x).^2;
crt = grt0./J;
crt(1,nt1+1:nt1+nt_inner)=crt(2,nt1+1:nt1+nt_inner);
%Due to singularity, the calculation of CRT at the origin was not successful. 
% Based on analytical approximation, its approximate constant was obtained. Therefore, it is added to facilitate the next step of differentiation
crt(index_x)=grt(index_x)./Jni(index_x)./R(index_x).^2;
ctt = gtt.*J./R.^2;

[dcrrdr,~]      =partial(crr,kdx,kdy,d1x,d1y,geometry);
[dcrtdr,dcrtdt] =partial(crt,kdx,kdy,d1x,d1y,geometry);
[~,dcttdt]      =partial(ctt,kdx,kdy,d1x,d1y,geometry);
%Seeking differentiation in another way to avoid the influence of X-point divergence on nearby derivative functions


[dJupR2dr,dJupR2dt]=partial(JupR2,kdx,kdy,d1x,d1y,geometry);
[dgrrdr,dgrrdt]=partial(grr,kdx,kdy,d1x,d1y,geometry);
[dgrtdr,dgrtdt]=partial(grt,kdx,kdy,d1x,d1y,geometry);
[dgttdr,dgttdt]=partial(gtt,kdx,kdy,d1x,d1y,geometry);
index_x=false(nr,nt);
range=15;
index_x(nr_down:nr,[t_min:nt1+1+range nt1+nt_inner+1-range:t_max])=true;%Take the reciprocal and derivative of the singular quantity in the lower area
%index_x(nr_inner+1-range:nr_inner+1+range,[nt1+1-range:nt1+1+range nt1+nt_inner+1-range:nt1+nt_inner+1+range])=true;
dcrrdr(index_x)=dgrrdr(index_x)./JupR2(index_x)-dJupR2dr(index_x).*grr(index_x)./JupR2(index_x).^2;
dcttdt(index_x)=dgttdt(index_x)./JupR2(index_x)-dJupR2dt(index_x).*gtt(index_x)./JupR2(index_x).^2;
dcrtdr(index_x)=dgrtdr(index_x)./JupR2(index_x)-dJupR2dr(index_x).*grt(index_x)./JupR2(index_x).^2;
dcrtdt(index_x)=dgrtdt(index_x)./JupR2(index_x)-dJupR2dt(index_x).*grt(index_x)./JupR2(index_x).^2;

M = struct( ...
    'dRdr',dRdr,'dZdr',dZdr ...
    ,'dRdt',dRdt,'dZdt',dZdt ...
    ,'J',J,'cJ',cJ,'Jni',Jni,'JupR2',JupR2...
    ,'crr',crr,'crt',crt,'ctt',ctt,'grr',grr,'grt',grt,'gtt',gtt  ...
    ,'dcrrdr',dcrrdr,'dcttdt',dcttdt ...
    ,'dcrtdr',dcrtdr,'dcrtdt',dcrtdt ...
    ,'dcrrdt',dcrrdr,'dcttdr',dcttdt);

end



function [Rr,Rt]=partial(R,kdx,kdy,d1x,d1y,geometry)
nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nt2=geometry.nt2;
nt=geometry.nt;
nr_inner=geometry.nr_inner;
nr_outer=geometry.nr_outer;
nr_down=geometry.nr_down;
nr=geometry.nr;
t_min=geometry.t_min;
t_max=geometry.t_max;

Rt=zeros(nr,nt);
Rr=zeros(nr,nt);
L=size(kdx,3);
for i=1:nr
    if i<nr_down
        for j=nt1+1:nt1+nt_inner
            for s=1:L
                Rr (i,j)=Rr (i,j)+d1x(i,j,s)*R(kdx(i,j,s),j);
                Rt (i,j)=Rt (i,j)+d1y(i,j,s)*R(i,kdy(i,j,s));
            end
        end
    else
        for j=t_min:t_max
            for s=1:L
                Rr (i,j)=Rr (i,j)+d1x(i,j,s)*R(kdx(i,j,s),j);
                Rt (i,j)=Rt (i,j)+d1y(i,j,s)*R(i,kdy(i,j,s));
            end
        end
    end
end
end

function Jnew=Eliminating_Singularities(J,err)
Jnew=J.*(abs(J)>=err)+err.*(abs(J)<err);
end

