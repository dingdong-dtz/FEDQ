function [f,A,B]=nabla2d_all_speed(rho,theta,L,Js,nd,bdy,geometry,R,Z,psi)

%************************************************************
% function [f,A,B,paras]=nabla2d_pol(x,y,L,Js,nd,bdy)
%   Set matrix A for linear equations  A * f = B
%   for 2D equations with peridic condition in y in [0,2*pi],
%   and a pole(origin) at x(1)=0,
%       i.e. a poloidal coordiantes with origin
%                      f = A\B;
%   constructed from differential equations:
%       L (f(x,y)) = Js(x,y),
%   where
%    L = L.dx2 * d_x^2 +L.dy2 * d_y^2 + L.dxdy *d_x d_y
%                + L.dx1 * d_x + L.dy1 * d_y + L.L0 ,
%   using mapping_uniform to get sequence number of points for derivatives
%    (i,j) ---> s   (unique)
%     1 <= i <= nx,  1 <= j <= (ny-1)
%      s  = 1,                if i==1   (only the first one unknown used)
%         = i +(j-1)*(nx-1),  else
%
%  Inputs:
%    x=x(1:nx)
%    y=x(1:ny),    y is in [0,2*pi]
%               y has to be uniform distributed and covering [0,2*pi]
%    ny,    grid number has to be 2^n+1,
%    L, structure for coefficients in the L operator
%  optional:
%    nd, order of accuracy 2^n, default =2
%     Js, for the right hand side of the equation
%   bdy, for boundary condition
%   bdy.psi_down
%
%  Ouputs:
%       f, f(1:nx,1:ny)
%          f is only solved when full information is ready
%     A,B,  for matrix form of equation A*f =B
%            A(lf,lf);   lf=(nx-1)*(ny-1)+1;
%            B(lf,1);
%    paras.
%      d2x,d1x,kdx,ldx,ixd,djxd
%      d2y,d1y,kdy,ldy
%            are for evaluate sencond order
%            and first order derivatives on the grids
%         using
%        d_x f |_j =  sum( dx1(j,:) .* f(kxd(j,:)) );
%        d_x^2 f |_j =  sum( dx2(j,:) .* f(kxd(j,:)) );
%************************************************************
nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nt2=geometry.nt2;
nt=geometry.nt;
nr_inner=geometry.nr_inner;
nr_outer=geometry.nr_outer;
nr_down=geometry.nr_down;
nr=geometry.nr;
t_min=geometry.t_min;
t_max=geometry.t_max;%
t_boundary=geometry.t_boundary;
lf=(nr_down-2)*nt_inner+1+(nr-nr_down+1)*(t_max-t_min+1);
[d2x,d1x,kdx,ldx,d2y,d1y,kdy,ldy]=getddiag_all(rho,theta,nd,geometry);
%----------------- set matrix
%----- set initial Matrix and default inputs
n_all=lf*((nd+1)^2+(nd+1)*2)+100;%稀疏矩阵预计具有的元素数
%下面以三个向量存储稀疏矩阵的存储位置和存贮值
Ai=zeros(n_all,1);
Aj=zeros(n_all,1);
Av=zeros(n_all,1);
B = zeros(lf,1);
S_ij=mapping_uniform(1:nr,1:nt,geometry);
%------ equation at middle core points
now=0;
for i=2:nr-1
    if i<=nr_down%第nr_down行外侧属于边界条件，需要在下方特殊处理
        jpre=nt1+1:nt1+nt_inner;
    else
        jpre=t_boundary(1,i)+1:t_boundary(2,i)-1;
    end
    for j=jpre
        if i==nr_inner+1&&(j==nt1+1||j==nt1+nt_inner+1)%x点需要在下方特殊处理
            continue
        end
        jkdy = kdy(i,j,:);
        jd2y = d2y(i,j,:);
        jd1y = d1y(i,j,:);
        
        ikdx = kdx(i,j,:);
        id2x = d2x(i,j,:);
        id1x = d1x(i,j,:);
        li=length(ikdx);
        lj=length(jkdy);
        %--------A
        s=S_ij(i,j);%  s row's of equations
        %------- Ax
        sx=S_ij(ikdx,j);
        Ai(now+1:now+li)=Ai(now+1:now+li)+s;
        Aj(now+1:now+li)=Aj(now+1:now+li)+sx(:);
        Av(now+1:now+li)=Av(now+1:now+li)+ L.dx2(i,j)*id2x(:) + L.dx1(i,j)*id1x(:);
        now=now+li;
        %------- Ay
        sy=S_ij(i,jkdy);
        Ai(now+1:now+lj)=Ai(now+1:now+lj)+s;
        Aj(now+1:now+lj)=Aj(now+1:now+lj)+sy(:);
        Av(now+1:now+lj)=Av(now+1:now+lj)+ L.dy2(i,j)*jd2y(:) + L.dy1(i,j)*jd1y(:);
        now=now+lj;
        %------- Axy
        sxy=S_ij(ikdx,jkdy);
        dxy=id1x(:)*jd1y(:)';
        lij=li*lj;
        Ai(now+1:now+lij)=Ai(now+1:now+lij)+s;
        Aj(now+1:now+lij)=Aj(now+1:now+lij)+sxy(:);
        Av(now+1:now+lij)=Av(now+1:now+lij)+ L.dxdy(i,j)*dxy(:);
        now=now+lij;
        %---------------- B
        B(s)= Js(i,j);
    end
end

%%
%-------------------- boundary condition and solution
% boundary position
bdy_out =S_ij(nr,t_min:t_max);
bdy_down=S_ij(nr_down,[t_min:nt1,nt1+nt_inner+1:t_max]);
l_out=length(bdy_out);
l_down=length(bdy_down);
Ai(now+1:now+l_out)=Ai(now+1:now+l_out)+bdy_out';
Aj(now+1:now+l_out)=Aj(now+1:now+l_out)+bdy_out';
Av(now+1:now+l_out)=Av(now+1:now+l_out)+1;
%B(bdy_out,1)=B(bdy_out,1)+bdy(nr);
B(bdy_out,1)=psi(nr,t_min:t_max)';
now=now+l_out;
Ai(now+1:now+l_down)=Ai(now+1:now+l_down)+bdy_down';
Aj(now+1:now+l_down)=Aj(now+1:now+l_down)+bdy_down';
Av(now+1:now+l_down)=Av(now+1:now+l_down)+1;
%B(bdy_down,1)=B(bdy_down,1)+bdy(nr_down);
B(bdy_down,1)=psi(nr_down,[t_min:nt1,nt1+nt_inner+1:t_max])';
now=now+l_down;
bdy_right=S_ij(nr_down+1:nr-1,t_min);
bdy_left =S_ij(nr_down+1:nr-1,t_max);
l_div=nr-nr_down-1;
Ai(now+1:now+l_div)=Ai(now+1:now+l_div)+bdy_right;
Aj(now+1:now+l_div)=Aj(now+1:now+l_div)+bdy_right;
Av(now+1:now+l_div)=Av(now+1:now+l_div)+1;
%B(bdy_right,1)=bdy(nr_down+1:nr-1);
B(bdy_right,1)=psi(nr_down+1:nr-1,t_min);
now=now+l_div;
Ai(now+1:now+l_div)=Ai(now+1:now+l_div)+bdy_left;
Aj(now+1:now+l_div)=Aj(now+1:now+l_div)+bdy_left;
Av(now+1:now+l_div)=Av(now+1:now+l_div)+1;
%B(bdy_left,1)=bdy(nr_down+1:nr-1);
B(bdy_left,1)=psi(nr_down+1:nr-1,t_max);
now=now+l_div;
%{
for k=1:length(bdy_out)%最外层磁面上的边界条件
    A(bdy_out(k),bdy_out(k))=1;
    B(bdy_out(k),1)=bdy(nr);
end
for k=1:length(bdy_down)%私有区磁面上的边界条件
    A(bdy_down(k),bdy_down(k))=1;
    B(bdy_down(k),1)=bdy(nr_down);
end

for i=nr_down+1:nr-1%偏滤器靶板处的边界条件
    s=S_ij(i,t_boundary(1:2,i)');
    A(s(1),s(1))=1;
    B(s(1),1)=bdy(i);
    A(s(2),s(2))=1;
    B(s(2),1)=bdy(i);
end
%}
A= sparse(Ai(1:now),Aj(1:now),Av(1:now),lf,lf,n_all);
%%
R0=R(1,nt1+1);
Z0=Z(1,nt1+1);
Rx=R(nr_inner+1,nt1+1);
Zx=Z(nr_inner+1,nt1+1);
[alpha,Length]=cart2pol(R0-Rx,Z0-Zx);
%-------------------- equation at axis
index_000=nt1+1;
[~,index_180]=min(abs(theta-pi));
[~,index_090]=min(abs(theta-pi/2));
[~,index_270]=min(abs(theta-pi/2*3));
index_sx_r=[3 2 1 2 3];index_sx_t=[index_090 index_090 index_090 index_270 index_270];
index_sy_r=[3 2 1 2 3];index_sy_t=[index_000 index_000 index_000 index_180 index_180];
%index_sx_r=[ 2 1 2 ];index_sx_t=[index_090 index_090 index_270];
%index_sy_r=[ 2 1 2 ];index_sy_t=[index_000 index_000 index_180];
x_axis=((R(sub2ind([nr,nt],index_sx_r,index_sx_t))-Rx)*sin(alpha)-(Z(sub2ind([nr,nt],index_sx_r,index_sx_t))-Zx)*cos(alpha))/Length;
y_axis=((R(sub2ind([nr,nt],index_sy_r,index_sy_t))-Rx)*cos(alpha)+(Z(sub2ind([nr,nt],index_sy_r,index_sy_t))-Zx)*sin(alpha))/Length;
x_d2x=fdcoeffF(2,0,x_axis);  %  d^2_x
x_d1x=fdcoeffF(1,0,x_axis);  %  d_x
x_d2y=fdcoeffF(2,0,y_axis);  %  d^2_x
x_d1y=fdcoeffF(1,0,y_axis);  %  d_x

sx=S_ij(sub2ind([nr,nt],index_sx_r,index_sx_t));%sub2ind是为了转化为线性索引
sy=S_ij(sub2ind([nr,nt],index_sy_r,index_sy_t));
s =S_ij(1,nt1+1);
A(s,sx) = A(s,sx) + x_d2x/Length^2 - x_d1x/Length/R0*sin(alpha);
A(s,sy) = A(s,sy) + x_d2y/Length^2 - x_d1y/Length/R0*cos(alpha);
B(s)=Js(1,nt1+1);

%%
%equation at x point
index_sx_r=[nr_inner+3 nr_inner+2 nr_inner+1 nr_inner+2 nr_inner+3];index_sx_t=[nt1+nt_inner+1 nt1+nt_inner+1 nt1+1 nt1+1 nt1+1];
index_sy_r=[nr_inner-1 nr_inner   nr_inner+1 nr_inner nr_inner-1  ];index_sy_t=[nt1+nt_inner+1 nt1+nt_inner+1 nt1+1 nt1+1 nt1+1];
%index_sx_r=[ nr_inner+2 nr_inner+1 nr_inner+2 ];index_sx_t=[ nt1+nt_inner+1 nt1+1 nt1+1];%三点格式
%index_sy_r=[ nr_inner   nr_inner+1 nr_inner   ];index_sy_t=[ nt1+nt_inner+1 nt1+1 nt1+1];
x_axis=((R(sub2ind([nr,nt],index_sx_r,index_sx_t))-Rx)*sin(alpha)-(Z(sub2ind([nr,nt],index_sx_r,index_sx_t))-Zx)*cos(alpha))/Length;
y_axis=((R(sub2ind([nr,nt],index_sy_r,index_sy_t))-Rx)*cos(alpha)+(Z(sub2ind([nr,nt],index_sy_r,index_sy_t))-Zx)*sin(alpha))/Length;
x_d2x=fdcoeffF(2,0,x_axis);  %  d^2_x
x_d1x=fdcoeffF(1,0,x_axis);  %  d_x
x_d2y=fdcoeffF(2,0,y_axis);  %  d^2_x
x_d1y=fdcoeffF(1,0,y_axis);  %  d_x

sx=S_ij(sub2ind([nr,nt],index_sx_r,index_sx_t));%sub2ind是为了转化为线性索引
sy=S_ij(sub2ind([nr,nt],index_sy_r,index_sy_t));
s =S_ij(nr_inner+1,nt1+1);
A(s,sx) = A(s,sx) + x_d2x/Length^2 - x_d1x/Length/Rx*sin(alpha);
A(s,sy) = A(s,sy) + x_d2y/Length^2 - x_d1y/Length/Rx*cos(alpha);
B(s)=Js(nr_inner+1,nt1+1);
%使矩阵中另一个也对应X点的位置与其一样
s_else=S_ij(nr_inner+1,nt1+nt_inner+1);
A(s_else,s)=-1;
A(s_else,s_else)=1;
B(s_else)=0;
%%
%---------- solusion, when all information of L(f)=Js is ready
f=zeros(nr,nt);
index_effect=false(nr,nt);
index_effect(1:nr_down-1,nt1+1:nt1+nt_inner)=true;
index_effect(nr_down:nr,t_min:t_max)=true;
f1=A\B;
f(index_effect)=f1(S_ij(index_effect));
%{
for i=1:nr
    if i<nr_down
        jpre=nt1+1:nt1+nt_inner;
    else
        jpre=t_boundary(1,i):t_boundary(2,i);
    end
    for j=jpre
        s=S_ij(i,j);
        f(i,j)=f1(s);
    end
end
%}
end