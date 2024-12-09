function [grr,drdR,drdZ,index_0,index_x,weight]=fit_axis_and_x(R,Z,rho,theta,geometry,X_n_fit)
%**************************************
% function [R0,Z0,npsi,psibm,psi_axis,sIp,r0,t0]=findaxis(eq,icheck)
% find magnetic axis
% 磁轴处等势面接近椭圆形，用内层磁面拟合确定磁轴位置
%*************************************
nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nr_down=geometry.nr_down;
nr_inner=geometry.nr_inner;
nr=geometry.nr;
nt=geometry.nt;
t_min=geometry.t_min;
t_max=geometry.t_max;

R0=R(1,nt1+1);
Z0=Z(1,nt1+1);
Rx=R(nr_inner+1,nt1+1);
Zx=Z(nr_inner+1,nt1+1);
[alpha,Length]=cart2pol(R0-Rx,Z0-Zx);
[Theta,Rho]=meshgrid(theta,rho);
npsi=Rho.^2;

drdR=zeros(nr,nt);
drdZ=zeros(nr,nt);
grr=zeros(nr,nt);

range_axis=5;
index_0=false(nr,nt);
index_0(2:range_axis,nt1+1:nt1+nt_inner)=true;
r02=2*(var(R(index_0))+var(Z(index_0)));%表示磁轴处拟合的取点范围
r0=sqrt(r02);
%{
figure
hold on
scatter(R(index_0),Z(index_0),8,'green','filled')
rectangle('Position',[R0-r0/sqrt(2),Z0-r0/sqrt(2),2*r0/sqrt(2),2*r0/sqrt(2)],'Curvature',[1 1]);
axis equal
%}
index_0=find((R-R0).^2+(Z-Z0).^2<r02*2);%重新以圆形区域选定参与拟合的网格点
weight_0=exp(-2*((R-R0).^2+(Z-Z0).^2)/r02);
[f_axis,gof] = fit([R(index_0),Z(index_0)],npsi(index_0),'poly22','Weights',weight_0(index_0));
%gof
theta_axis=(Theta(index_0)-pi+alpha);
dpsidr=f_axis.p20*cos(theta_axis).^2+f_axis.p02*sin(theta_axis).^2+f_axis.p11*sin(theta_axis).*cos(theta_axis);
drdR(index_0)=(f_axis.p20*cos(theta_axis)+0.5*f_axis.p11*sin(theta_axis))./sqrt(dpsidr);
drdZ(index_0)=(f_axis.p02*sin(theta_axis)+0.5*f_axis.p11*cos(theta_axis))./sqrt(dpsidr);
grr(index_0)=drdR(index_0).^2+drdZ(index_0).^2;

index_x=false(nr,nt);
index_x(nr_inner+1-X_n_fit:nr_inner+1+X_n_fit,[nt1+1-X_n_fit:nt1+1+X_n_fit nt1+nt_inner+1-X_n_fit:nt1+nt_inner+1+X_n_fit])=true;
near_R=R(index_x);  near_Z=Z(index_x); 
rx2=(var(near_R)+var(near_Z));%表示拟合区域的平均半径
rx=sqrt(rx2);
%{
figure
hold on
scatter(near_R,near_Z)
rectangle('Position',[Rx-rx/sqrt(2),Zx-rx/sqrt(2),2*rx/sqrt(2),2*rx/sqrt(2)],'Curvature',[1 1]);
axis equal
%}
index_x=find((R-Rx).^2+(Z-Zx).^2<2*rx2);%重新以圆形区域选定参与拟合的网格点
weight_x=exp(-2*((R-Rx).^2+(Z-Zx).^2)/rx2);
[f_x,~] = fit([R(index_x),Z(index_x)],npsi(index_x),'poly22','Weights',weight_x(index_x));
%gof
drdR(index_x)=(f_x.p20*R(index_x)+0.5*f_x.p11*Z(index_x)+0.5*f_x.p10)./Rho(index_x);
drdZ(index_x)=(f_x.p02*Z(index_x)+0.5*f_x.p11*R(index_x)+0.5*f_x.p01)./Rho(index_x);
grr(index_x)=drdR(index_x).^2+drdZ(index_x).^2;
weight_x=exp(-2*((R-R0).^2+(Z-Z0).^2)/r02)+exp(-2*(((R-Rx).^2+(Z-Zx).^2)/rx2).^2);

weight=weight_0+weight_x;%全空间以待用拟合更新的权重
end