function [Rx,Zx,npsi_x,npsi_fit,SV,M,T,r] = fit_xpoint_new(Rx_0,Zx_0,R,Z,npsi,geometry,icheck,X_n_fit)
%[Rx,Zx,npsi_x,npsi_fit,SV,M,T,index_x] = fit_xpoint_new(Rx_0,Zx_0,R,Z,npsi,geometry,icheck,X_n_range)
%Rx_0,Zx_0  上一步程序中的X点位置
%对x点做拟合，找到x点位置和磁通值
%X_n_range  表示在拟合时选取十字线两侧各X_n_range条线相交处的网格点参与拟合
%index_x    是x点附近做拟合时选用的点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nt=geometry.nt; nt1=geometry.nt1;   nt_inner=geometry.nt_inner;
nr=geometry.nr; nr_inner=geometry.nr_inner;
index_x=false(nr,nt);
index_x(nr_inner+1-X_n_fit:nr_inner+1+X_n_fit,[nt1+1-X_n_fit:nt1+1+X_n_fit nt1+nt_inner+1-X_n_fit:nt1+nt_inner+1+X_n_fit])=true;
near_R=R(index_x);  near_Z=Z(index_x);  near_npsi=npsi(index_x);
r2=var(near_R)+var(near_Z);%表示拟合区域的平均半径
r=sqrt(r2);
index=find((R-Rx_0).^2+(Z-Zx_0).^2<r2);
near_R=R(index);
near_Z=Z(index);
near_npsi=npsi(index);
%{
figure
hold on
scatter(near_R,near_Z)
rectangle('Position',[Rx_0-r,Zx_0-r,2*r,2*r],'Curvature',[1 1]);
rectangle('Position',[Rx_0-r/1.4,Zx_0-r/1.4,2*r/1.4,2*r/1.4],'Curvature',[1 1]);
axis equal
%}
%{
[rho_n,theta_n]=meshgrid(1:nt,1:nr);
rho_n=abs(rho_n-(nr_inner+1));
theta_n=min(abs(theta_n-(nt1+1)),abs(theta_n-(nt1+1+nt_inner)));
weight=exp(-(rho_n(index_x).^2+theta_n(index_x).^2)/X_n_fit^2);
%}
weight=exp(-2*((near_R-Rx_0).^2+(near_Z-Zx_0).^2)/r2);%使靠近X点区域误差更小，因为越靠近X点越符合二元二次多项式
%[sf,~] = fit([near_R, near_Z],near_npsi,'poly22');
[sf,gof] = fit([near_R, near_Z],near_npsi,'poly22','Weights',weight);
%gof
%下面对于拟合出的结果做配方
M=[sf.p20,sf.p11/2;
    sf.p11/2,sf.p02];
T=[sf.p10 sf.p01];
[V,D]=eig(M);%M*V = V*D
H=V;
TH=T*H;

npsi_x=sf.p00-TH(1)^2/4/D(1,1)-TH(2)^2/4/D(2,2);
X_Point=H*[-TH(1)/2/D(1,1);-TH(2)/2/D(2,2)];
Rx=X_Point(1);
Zx=X_Point(2);
S=[sqrt(abs(D(1,1))),sqrt(abs(D(2,2)));sqrt(abs(D(1,1))),-sqrt(abs(D(2,2)))]/(abs(D(1,1))+abs(D(2,2)));
SV=S*V';

npsi_fit_diff=(sf(R,Z)-npsi).*exp(-((R-Rx).^2+(Z-Zx).^2)/(r2/2));

disp(['最大拟合偏差',num2str(max_geometry(abs(npsi_fit_diff),geometry),'%10.1e')]);
npsi_fit=npsi+npsi_fit_diff;

if icheck
    plot_matrix_delta(R,Z,npsi_fit_diff,'npsi_fit_diff',geometry,[])
    npsi_pri=npsi(geometry.nr_down,geometry.nt1);
    npsi_lcs=npsi(end,geometry.nt1);
    plot_matrix(R,Z,npsi,'psi_fit',geometry,[],linspace(npsi_pri,npsi_lcs,25),0)
    plot_matrix(R,Z,npsi_fit,[],geometry,[],linspace(npsi_pri,npsi_lcs,25),0)
end
end