function [gtt,t_R,t_Z]  = gtt_theory(R,Z,theta,rho,geometry)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nt2=geometry.nt2;
nt=geometry.nt;
nr_inner=geometry.nr_inner;
nr_outer=geometry.nr_outer;
nr_down=geometry.nr_down;
nr=geometry.nr;
t_boundary=geometry.t_boundary;
lf=(nr_down-2)*nt_inner+1+sum(t_boundary(3,nr_down:nr));% number of unkonwns
t_min=geometry.t_min;
t_max=geometry.t_max;
[Theta,Rho]=meshgrid(theta,rho);
R0=R(1,nt1+1);
Z0=Z(1,nt1+1);
Rx=R(nr_inner+1,nt1+1);
Zx=Z(nr_inner+1,nt1+1);
[alpha,Length]=cart2pol(R0-Rx,Z0-Zx);
[x,y] = rotate(R,Z,R0,Z0,alpha,Length);

gtt=zeros(nr,nt);
t_y=zeros(nr,nt);
t_x=zeros(nr,nt);
ifdown=zeros(size(y));
ifdown(nr_down:nr,[t_min:nt1 nt1+nt_inner+1:t_max])=1;
xbar=pi/2*x;
ybar=pi/2*(y+2*ifdown);

index1=false(nr,nt);
index1(1:nr,nt1+1:nt1+nt_inner)=true;
index1=index1&(Theta<pi/2|Theta>pi/2*3);%中部磁轴处

gtt(index1)=pi^2/4/Length^2*(xbar(index1).^2.*cos(Theta(index1)).^4./sin(ybar(index1)).^4+sin(Theta(index1)).^2.*cos(Theta(index1)).^2./xbar(index1).^2);
t_y(index1)=pi/2/Length*(xbar(index1).*cos(Theta(index1)).^2./sin(ybar(index1)).^2);
t_x(index1)=pi/2/Length*(sin(Theta(index1)).*cos(Theta(index1))./xbar(index1));
%%
index1=false(nr,nt);
index1(1:nr,nt1+1:nt1+nt_inner)=true;
index1=index1&(Theta<pi/2&Theta>pi/2-pi/9|Theta>pi/2*3&Theta<pi/2*3+pi/9);%中部y=0附近
gtt(index1)=pi^2/4/Length^2*(xbar(index1).^(-2).*sin(Theta(index1)).^4./cos(ybar(index1)).^4+...
    sin(Theta(index1)).^4./xbar(index1).^4.*tan(ybar(index1)).^2);
t_y(index1)=pi/2/Length*(xbar(index1).^(-1).*sin(Theta(index1)).^2./cos(ybar(index1)).^2);
t_x(index1)=pi/2/Length*(-sin(Theta(index1)).^2./xbar(index1).^2.*tan(ybar(index1)));
%%
index1=false(nr,nt);
index1(nr_down:nr,[t_min:nt1 nt1+nt_inner+1:t_max])=true;%下方
gtt(index1)=pi^2/4/Length^2*(xbar(index1).^2.*cos(Theta(index1)).^4./sin(ybar(index1)).^4+sin(Theta(index1)).^2.*cos(Theta(index1)).^2./xbar(index1).^2);
t_y(index1)=pi/2/Length*(xbar(index1).*cos(Theta(index1)).^2./sin(ybar(index1)).^2);
t_x(index1)=pi/2/Length*(sin(Theta(index1)).*cos(Theta(index1))./xbar(index1));
%%

index1=false(nr,nt);
index1(1:nr_inner,nt1+1)=true;index1(nr_down:nr_inner,nt1+nt_inner+1)=true;
gtt(index1)=pi^2/4/Length^2*(cos(Theta(index1)).^2.*sin(Theta(index1)).^2./sin(ybar(index1)).^2./cos(ybar(index1)).^2+...
    cos(Theta(index1)).^4./tan(ybar(index1)).^2);
t_y(index1)=pi/2/Length*(-cos(Theta(index1)).*sin(Theta(index1))./sin(ybar(index1))./cos(ybar(index1)));
t_x(index1)=pi/2/Length*(-cos(Theta(index1)).^2./tan(ybar(index1)));
%%
gtt(nr_inner+1,[nt1+1,nt1+nt_inner+1])=0;%X点
t_y(nr_inner+1,[nt1+1,nt1+nt_inner+1])=0;
t_x(nr_inner+1,[nt1+1,nt1+nt_inner+1])=0;
%%
index1=false(nr,nt);
index1(1:nr,nt1+1:nt1+nt_inner)=true;
index1=index1&(Theta>=pi/2&Theta<=pi/2*3);%上方区域
index1(1:10,nt1+1:nt1+nt_inner)=true;
gtt(index1)=1/Length^2./(x(index1).^2+y(index1).^2);
t_y(index1)=1/Length*(x(index1))./(x(index1).^2+y(index1).^2);
t_x(index1)=1/Length*(-y(index1))./(x(index1).^2+y(index1).^2);


t_R=+sin(alpha)*t_x+cos(alpha)*t_y;
t_Z=-cos(alpha)*t_x+sin(alpha)*t_y;

%delta=t_y.^2+t_x.^2-gtt;

end


function [x,y] = rotate(R,Z,R0,Z0,alpha,L)
%[x,y] = rotate(R,Z,R0,Z0,alpha,L)
%将坐标系顺时针旋转pi/2-α,以（R0,Z0）为中心并缩小L倍
x=((R-R0)*sin(alpha)-(Z-Z0)*cos(alpha))/L;
y=((R-R0)*cos(alpha)+(Z-Z0)*sin(alpha))/L;
end
function [R,Z] = rotate_inverse(x,y,R0,Z0,alpha,L)
%[R,Z] = rotate_inverse(x,y,R0,Z0,alpha,L)
%以（0,0）为中心放大L倍,再将坐标系逆时针旋转pi/2-α
R=(cos(alpha)*y+sin(alpha)*x)*L+R0;
Z=(sin(alpha)*y-cos(alpha)*x)*L+Z0;
end
