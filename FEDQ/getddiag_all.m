function [d2x,d1x,kdx,ldx,d2y,d1y,kdy,ldy]=getddiag_all(r,t,nd,geometry)
%*****************************************************
%  get diagnal parts of the matrix form of differential operator,
%  it has
%      d_x f |_j =  sum( d1(j,:) .* f(kd(j,:)) );
%    d_x^2 f |_j =  sum( d2(j,:) .* f(kd(j,:)) );
% Outputs:
%     d2, the coefficient for the d_x^2
%     d1, the coefficient for the d_x
%     kd, = j + ld, the sequence number of the grids for field evaluation
% Inputs:
%     r(:), rho grids;  t(:),theta grids
%     nd, order of accuracy 2^n
%     groemtry：网格区域划分与有效格点边界
%*****************************************************
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
d2x=zeros(nr,nt,nd+1);
d1x=zeros(nr,nt,nd+1);
ldx=zeros(nr,nt,nd+1);
kdx=zeros(nr,nt,nd+1);
d2y=zeros(nr,nt,nd+1);
d1y=zeros(nr,nt,nd+1);
ldy=zeros(nr,nt,nd+1);
kdy=zeros(nr,nt,nd+1);

%----------------- set matrix
[d2x_up,d1x_up,kdx_up,ldx_up]=getddiag(r,nd,1,0);%上方rho方向的导数系数
[d2x0,d1x0,kdx0,ldx0]=getddiag(r(nr_down:nr),nd,1,0);%下方rho方向的导数系数
d2x_down=d2x0;                          d1x_down=d1x0;
kdx_down=kdx0+nr_down-1;                ldx_down=ldx0;%极向坐标都是从nt1+1开始，只包含nt_inner个点
%将矩阵转化为正确的格式方便之后在循环中赋值
[d2y0,d1y0,kdy0,ldy0]=getddiag(t(nt1+1:nt1+nt_inner),nd,2,1);
d2y_inner=d2y0;                          d1y_inner=d1y0;
kdy_inner=kdy0+nt1;                      ldy_inner=ldy0;%极向坐标都是从nt1+1开始，只包含nt_inner个点
[d2y0,d1y0,kdy0,ldy0]=getddiag([t(t_min:nt1),t(nt1+nt_inner+1:t_max)-2*pi],nd,2,0);
d2y_pri=d2y0;                                        d1y_pri=d1y0;
kdy_pri=kdy0+t_min-1+(kdy0>nt1-t_min+1)*nt_inner;    ldy_pri=ldy0;%极向坐标都是从t_min开始，只包含nt_inner个点
[d2y0,d1y0,kdy0,ldy0]=getddiag(t(t_min:t_max),nd,2,0);
d2y_SOL=d2y0;                            d1y_SOL=d1y0;
kdy_SOL=kdy0+t_min-1;                    ldy_SOL=ldy0;%极向坐标都是从t_min开始，只包含nt_inner个点

for i=2:nr_inner
    kk=nt1+1:nt1+nt_inner;%本部分表示的是核心区域等离子体
    ldy(i,kk,:) = ldy_inner;
    kdy(i,kk,:) = kdy_inner;
    d2y(i,kk,:) = d2y_inner;
    d1y(i,kk,:) = d1y_inner;
    if i>=nr_down%-- extend t before evaluation 对于私有区域的点
        kk=[t_min:nt1,nt1+nt_inner+1:t_max];
        d2y(i,kk,:)=d2y_pri;     d1y(i,kk,:)=d1y_pri;
        kdy(i,kk,:)=kdy_pri;     ldy(i,kk,:)=ldy_pri;
    end
end
for i=nr_inner+1:(nr-1)
    kk=(t_min:t_max);
    d2y(i,kk,:)=d2y_SOL;     d1y(i,kk,:)=d1y_SOL;
    kdy(i,kk,:)=kdy_SOL;     ldy(i,kk,:)=ldy_SOL;
end
for j=nt1+1:nt1+nt_inner%本部分表示的是x点上方区域
    d1x(:,j,:) = d1x_up;     d2x(:,j,:) = d2x_up;
    kdx(:,j,:) = kdx_up;     ldx(:,j,:) = ldx_up;
end
for j=[1:nt1,nt1+nt_inner+1:nt]%本部分表示的是x点下方区域
    kk=nr_down:nr;
    d2x(kk,j,:)=d2x_down;     d1x(kk,j,:)=d1x_down;
    kdx(kk,j,:)=kdx_down;     ldx(kk,j,:)=ldx_down;
end


%下面部分将分界线和十字架的X点处的差分格式更新一下
j=nt1+1;
[d2x_x,d1x_x,kdx_x,ldx_x]=getddiag(r(1:nr_inner+1),nd,1,0);%上方rho方向的导数系数
d2x(1:nr_inner+1,j,:)=d2x_x;
d1x(1:nr_inner+1,j,:)=d1x_x;
kdx(1:nr_inner+1,j,:)=kdx_x;
ldx(1:nr_inner+1,j,:)=ldx_x;
[d2x_x,d1x_x,kdx_x,ldx_x]=getddiag(r(nr_inner+1:nr),nd,1,0);%上方rho方向的导数系数
d2x(nr_inner+1:nr,j,:)=d2x_x;
d1x(nr_inner+1:nr,j,:)=d1x_x;
kdx(nr_inner+1:nr,j,:)=kdx_x+nr_inner;
ldx(nr_inner+1:nr,j,:)=ldx_x;
j=nt1+nt_inner+1;
[d2x_x,d1x_x,kdx_x,ldx_x]=getddiag(r(nr_down:nr_inner+1),nd,1,0);%上方rho方向的导数系数
d2x(nr_down:nr_inner+1,j,:)=d2x_x;
d1x(nr_down:nr_inner+1,j,:)=d1x_x;
kdx(nr_down:nr_inner+1,j,:)=kdx_x+nr_down-1;
ldx(nr_down:nr_inner+1,j,:)=ldx_x;
[d2x_x,d1x_x,kdx_x,ldx_x]=getddiag(r(nr_inner+1:nr),nd,1,0);%上方rho方向的导数系数
d2x(nr_inner+1:nr,j,:)=d2x_x;
d1x(nr_inner+1:nr,j,:)=d1x_x;
kdx(nr_inner+1:nr,j,:)=kdx_x+nr_inner;
ldx(nr_inner+1:nr,j,:)=ldx_x;
i=nr_inner+1;
[d2y0,d1y0,kdy0,ldy0]=getddiag(t(nt1+1:nt1+nt_inner),nd,2,0);%分界面在X点处不光滑
d2y(i,nt1+1:nt1+nt_inner,:)=d2y0;                          d1y(i,nt1+1:nt1+nt_inner,:)=d1y0;
kdy(i,nt1+1:nt1+nt_inner,:)=kdy0+nt1;                      ldy(i,nt1+1:nt1+nt_inner,:)=ldy0;%极向坐标都是从nt1+1开始，只包含nt_inner个点
[d2y0,d1y0,kdy0,ldy0]=getddiag(t(t_min:nt1),nd,2,0);
d2y(i,t_min:nt1,:)=d2y0;                          d1y(i,t_min:nt1,:)=d1y0;
kdy(i,t_min:nt1,:)=kdy0+t_min-1;                  ldy(i,t_min:nt1,:)=ldy0;%极向坐标都是从nt1+1开始，只包含nt_inner个点
[d2y0,d1y0,kdy0,ldy0]=getddiag(t(nt1+nt_inner+1:t_max)-2*pi,nd,2,0);
d2y(i,nt1+nt_inner+1:t_max,:)=d2y0;                 d1y(i,nt1+nt_inner+1:t_max,:)=d1y0;
kdy(i,nt1+nt_inner+1:t_max,:)=kdy0+nt1+nt_inner;    ldy(i,nt1+nt_inner+1:t_max,:)=ldy0;
end
function [d2,d1,kd,ld]=getddiag(x,nd,x_or_y,ifperiodic)
%*****************************************************
%  get diagnal parts of the matrix form of
%     differential operator,
%     it has
%      d_x f |_j =  sum( d1(j,:) .* f(kd(j,:)) );
%      d_x^2 f |_j =  sum( d2(j,:) .* f(kd(j,:)) );
% Outputs:
%     d2, the coefficient for the d_x^2
%     d1, the coefficient for the d_x
%     kd, = j + ld, the sequence number of the grids
%         for field evaluation
% Inputs:
%     x(:), x grids
%     nd, order of accuracy 2^n
%     x_or_y:确定此次调用求解的是ρ方向还是θ方向
%     ifperiodic：确定此次调用求解的模块是否为周期性边界条件，=0时表示非周期条件
%*****************************************************
nx=length(x);
ng=floor(nd/2);
kd0=-ng +(0:nd);
d2=zeros(nx,nd+1);
d1=zeros(nx,nd+1);
ld=zeros(nx,nd+1);
kd=zeros(nx,nd+1);
if ifperiodic==0%非周期条件
    for j=1:nx
        jkd = snearby(j,kd0,nx);%寻找j位置附近的几个点
        [d2(j,:),d1(j,:)]=getd(x(j),x(jkd));%求出x（j）处的导数
        kd(j,:) = jkd;%记录求导数用的哪些点
        ld(j,:) = jkd-j;%差分格式中使用的点和原先的点的关系
    end
else
    for j=1:nx
        jkd=mod(kd0+j-1,nx)+1;
        [d2(j,:),d1(j,:)]=getd(x(j),x(jkd)+2*pi*floor((kd0+j-1)/nx));%求出x（j）处的导数
        kd(j,:) = jkd;%记录求导数用的哪些点
        ld(j,:) = jkd-j;%差分格式中使用的点和原先的点的关系
    end
end
if x_or_y==1%求r方向的差分时
    d2=permute(d2,[1,3,2]);d1=permute(d1,[1,3,2]);kd=permute(kd,[1,3,2]);ld=permute(ld,[1,3,2]);
else%求t方向的差分时
    d2=permute(d2,[3,1,2]);d1=permute(d1,[3,1,2]);kd=permute(kd,[3,1,2]);ld=permute(ld,[3,1,2]);
end
end

function [d2,d1]=getd(xbar,x)
%*****************************************************
%  get coefficients for 1d differential operators
%    xbar is in the range of x
%    and the length of x determine the accuracy
%    for xbar=x_j,   x=x(j-1:j+1),
%         it is second order accuracy
%*****************************************************
d2=fdcoeffF(2,xbar,x);  %  d^2_x
d1=fdcoeffF(1,xbar,x);  %  d_x
end
function jl = snearby(j,l0,nx)
%****************************************************
%   jl= j + l0
%    shift grids for ending points
%****************************************************
jl=j+l0;
%--- for ending points
if jl(1)<=0
    jl=jl-jl(1)+1;
elseif jl(end)>=nx
    jl=jl-jl(end)+nx;
end
end
function c = fdcoeffF(k,xbar,x)

n = length(x);
if k >= n
    error('*** length(x) must be larger than k')
end

m = k;   % change to m=n-1 if you want to compute coefficients for all
% possible derivatives.  Then modify to output all of C.
c1 = 1;
c4 = x(1) - xbar;
C = zeros(n-1,m+1);
C(1,1) = 1;
for i=1:n-1
    i1 = i+1;
    mn = min(i,m);
    c2 = 1;
    c5 = c4;
    c4 = x(i1) - xbar;
    for j=0:i-1
        j1 = j+1;
        c3 = x(i1) - x(j1);
        c2 = c2*c3;
        if j==i-1
            for s=mn:-1:1
                s1 = s+1;
                C(i1,s1) = c1*(s*C(i1-1,s1-1) - c5*C(i1-1,s1))/c2;
            end
            C(i1,1) = -c1*c5*C(i1-1,1)/c2;
        end
        for s=mn:-1:1
            s1 = s+1;
            C(j1,s1) = (c4*C(j1,s1) - s*C(j1,s1-1))/c3;
        end
        C(j1,1) = c4*C(j1,1)/c3;
    end
    c1 = c2;
end

c = C(:,end)';            % last column of c gives desired row vector
end