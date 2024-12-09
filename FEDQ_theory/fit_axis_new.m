function [R0,Z0,psi_axis,M,npsi_fit,psib]=fit_axis_new(R,Z,psi,geometry,rho)
%**************************************
% function [R0,Z0,npsi,psibm,psi_axis,sIp,r0,t0]=findaxis(eq,icheck)
% find magnetic axis
% 磁轴处等势面接近椭圆形，用内层磁面拟合确定磁轴位置
% 具体的选定内层五个磁面做拟合
%*************************************
nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nr_down=geometry.nr_down;
nr=geometry.nr;
nt=geometry.nt;
t_min=geometry.t_min;
t_max=geometry.t_max;

psi_lcs=psi(end,nt1);
psi_axis0=-max_geometry(-psi,geometry);
R0=R(1,nt1+1);
Z0=Z(1,nt1+1);
psi_new=psi(2:nr_down,nt1+1:nt1+nt_inner);  psi_new=[psi(1,nt1+1);psi_new(:)];
R_new=R(2:nr_down,nt1+1:nt1+nt_inner);      R_new=[R(1,nt1+1); R_new(:)];
Z_new=Z(2:nr_down,nt1+1:nt1+nt_inner);      Z_new=[Z(1,nt1+1); Z_new(:)];%只取核心区域的点来计算
r2_for_0=(R_new(5,1)-R0)^2+(Z_new(5,1)-Z0)^2;
psi_pre=(psi_lcs-psi_axis0)*rho(6)^2+psi_axis0;%因为前五层采用此方式更新，故就选取涉及前五层的值
index=find(psi_new<psi_pre);
ssss=0;
while length(index)<nt_inner*3
    ssss=ssss+1;
    index=find(psi_new<rho(5+ssss)^2*(psi_lcs-psi_axis0)+psi_axis0);
end
sf = fit([R_new(index),Z_new(index)],psi_new(index),'poly22');
M=[sf.p20,sf.p11/2;
    sf.p11/2,sf.p02];
T=[sf.p10 sf.p01];
[V,D]=eig(M);%M*V = V*D
%e1=V(:,1)/(V(:,1)'*V(:,1));%做正交归一化，但是函数生成的V已经是正交归一的，这样做没有价值
%e2=V(:,2)/(V(:,2)'*V(:,2));
%H=[e1,e2];
H=V;
TH=T*H;
psi_axis=sf.p00-TH(1)^2/4/D(1,1)-TH(2)^2/4/D(2,2);
X_Point=H*[-TH(1)/2/D(1,1);-TH(2)/2/D(2,2)];
R0=X_Point(1);
Z0=X_Point(2);

% normalize polidal flux
psib= psi_lcs-psi_axis;
psi_fit_diff=(sf(R,Z)-psi).*exp(-((R-R0).^2+(Z-Z0).^2)/(r2_for_0));
psi_fit=psi+psi_fit_diff;
npsi_fit=max((psi_fit-psi_axis)./psib,0); % normalized polidal flux,且使得最小值不小于0，以免开根号出现虚数
M=M/psib;
end