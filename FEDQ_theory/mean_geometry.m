function Mmean = mean_geometry(M,geometry,J,R,Z,rho,theta)
%   Mmean = mean_geometry(M,geometry,eq)
%   该函数用于求解整个等离子体空间的平均值，注意是体空间的平均值，也就是考虑了大环半径的效应
nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nt=geometry.nt;
nr_inner=geometry.nr_inner;
nr_down=geometry.nr_down;
nr=geometry.nr;
t_min=geometry.t_min;
t_max=geometry.t_max;
%筛选出有效数据区域
index=false(nr,nt);
index(1:nr_down-1,nt1+1:nt1+nt_inner)=true;
index(nr_down:nr,t_min:t_max)=true;


%求出\rho\theta空间中面积元
delta_rho=diff(rho);
delta_rho=1/2*([delta_rho,0]+[0,delta_rho]);
theta_inner=theta(nt1+1:nt1+nt_inner);
delta_theta_inner=diff(theta_inner);
delta_theta_inner=([delta_theta_inner,delta_theta_inner(1)]+[delta_theta_inner(end),delta_theta_inner])/2;
theta_SOL=theta(t_min:t_max);
delta_theta_SOL=diff(theta_SOL);
delta_theta_SOL=1/2*([delta_theta_SOL,0]+[0,delta_theta_SOL]);
theta_pri=[theta(t_min:nt1),theta(nt1+nt_inner+1:t_max)-2*pi];
delta_theta_pri=diff(theta_pri);
delta_theta_pri=1/2*([delta_theta_pri,0]+[0,delta_theta_pri]);
delta_S_rt=zeros(nr,nt);
delta_S_rt(1:nr_inner,nt1+1:nt1+nt_inner)=delta_rho(1:nr_inner)'*delta_theta_inner;
delta_S_rt(nr_inner+1:nr,t_min:t_max)=delta_rho(nr_inner+1:nr)'*delta_theta_SOL;
delta_S_rt(nr_inner+1:nr,[t_min:nt1,nt1+nt_inner+1:t_max])=delta_rho(nr_inner+1:nr)'*delta_theta_pri;
%转化为实际体空间的体积元
delta_S_RZ=delta_S_rt.*J;
%求出X点处对应的小面元的面积
rho_near_x=1+[nr_inner-1,nr_inner-1,nr_inner-1,nr_inner,nr_inner+1,nr_inner+1,nr_inner+1,nr_inner,nr_inner-1,nr_inner-1,nr_inner-1,nr_inner,nr_inner+1,nr_inner+1,nr_inner+1,nr_inner,nr_inner-1];
nl=nt1+nt_inner;
theta_near_x=[nt1+2,nt1+1,nl,nl,nl,nl+1,nl+2,nl+2,nl+2,nl+1,nt1,nt1,nt1,nt1+1,nt1+2,nt1+2,nt1+2];
R_near=R(sub2ind(size(R),rho_near_x,theta_near_x));
Z_near=Z(sub2ind(size(Z),rho_near_x,theta_near_x));
Sx=polyarea(R_near,Z_near)/8;
delta_S_RZ(nr_inner+1,nt1+1)=Sx;
delta_S_RZ(nr_inner+1,nt1+nt_inner+1)=Sx;
%求出平均值
Mmean=sum(M(index).*delta_S_RZ(index))/sum(index(index).*delta_S_RZ(index));
end