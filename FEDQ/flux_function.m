function [Pressure,g_RBt,Bt]=flux_function(gpsi,Pprime,FFprime,gpsi_pri,Pprime_pri,FFprime_pri,R,psi,geometry,psi_x)
%*****************************************************
% function Js=evalJs(gpsi,Pprime,FFprime,R,nrho,npsi)
%  evaluate right hand side of GS equaiton
%       Js = - [ R^2 * mu0 * P' + FF'];
% inputs:
%       gpsi,       normalized psi_p grids for P' and FF'
%       Pprime,     P'
%       FFprime     FF'
%       gpsi_pri
%       Pprime_pri
%       FFprime_pri
%       R,          major radius of the surfaces(按矩阵位置排列的)
%       psi,        psi_p (按矩阵位置排列的)psi要与gpsi对应，要么都是归一化的，要么都是不归一化的
%******************************************************

if ~exist('Pedg','var')
    Pedg=0;
end
if ~exist('gedg','var')
    gedg=4;%表示磁轴处磁场强度
end
P0=cumtrapz(gpsi,Pprime);
P=P0-P0(end)+Pedg;
Ppri0=cumtrapz(gpsi_pri,Pprime_pri);%必须保证gpsi_pri也是从
Ppri=Ppri0-P0(end)+Pedg;

g0=cumtrapz(gpsi,FFprime);
g=g0-g0(1)+gedg;
gpri0=cumtrapz(gpsi_pri,FFprime_pri);%必须保证gpsi_pri也是从
gpri=gpri0-gpri0(1)+gedg;

P_x=interp1(gpsi,P,psi_x,'pchip');
P_x0=interp1(gpsi,Ppri,psi_x,'pchip');
Ppri=Ppri-P_x0+P_x;

g_x=interp1(gpsi,g,psi_x,'pchip');
g_x0=interp1(gpsi,gpri,psi_x,'pchip');
gpri=gpri-g_x0+g_x;


nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nt=geometry.nt;
nr_inner=geometry.nr_inner;
nr=geometry.nr;
t_min=geometry.t_min;
t_max=geometry.t_max;
Pressure=zeros(nr,nt);
g_RBt=zeros(nr,nt);
Pressure(1:nr_inner,nt1+1:nt1+nt_inner)=interp1(gpsi,P,psi(1:nr_inner,nt1+1:nt1+nt_inner),'pchip');  
g_RBt(1:nr_inner,nt1+1:nt1+nt_inner)=interp1(gpsi,g,psi(1:nr_inner,nt1+1:nt1+nt_inner),'pchip');
Pressure(nr_inner+1:end,:)=interp1(gpsi,P,psi(nr_inner+1:end,:),'pchip'); 
g_RBt(nr_inner+1:end,:)=interp1(gpsi,g,psi(nr_inner+1:end,:),'pchip');
Pressure(1:nr_inner,[t_min:nt1 nt1+nt_inner+1:t_max])=interp1(gpsi_pri,Ppri,psi(1:nr_inner,[t_min:nt1 nt1+nt_inner+1:t_max]),'pchip');  % Pprime is defined in npsi
g_RBt(1:nr_inner,[t_min:nt1 nt1+nt_inner+1:t_max])=interp1(gpsi_pri,gpri,psi(1:nr_inner,[t_min:nt1 nt1+nt_inner+1:t_max]),'pchip');
Bt=g_RBt./R;
end
