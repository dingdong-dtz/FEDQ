function [d2x,d1x,kdx,d2y,d1y,kdy]=diff_me_steady_change(rho,theta,number,order,geometry)
%*****************************************************
%  get diagnal parts of the matrix form of differential operator,
%  it has
%      d_x f |_j =  sum( d1(j,:) .* f(kd(j,:)) );
%    d_x^2 f |_j =  sum( d2(j,:) .* f(kd(j,:)) );
% Outputs:
%     d2, 二阶导的差分格式系数
%     d1, 一阶导的差分格式系数
%     kd, 差分格式中使用的格点的位置
% Inputs:
%     rho(:), rho grids;  theta(:),theta grids
%     number, 求每个差分导数时使用的格点数目
%     groemtry  网格区域划分与有效格点边界
%     order     求导计算中的阶数
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
d2x=zeros(nr,nt,number);
d1x=zeros(nr,nt,number);
kdx=zeros(nr,nt,number);
d2y=zeros(nr,nt,number);
d1y=zeros(nr,nt,number);
kdy=zeros(nr,nt,number);
%----------------- set matrix
[d2x_up,d1x_up,kdx_up]=getddiag_steady(rho,number,order,1,0);%上方rho方向的导数系数
[d2x0,d1x0,kdx0]=getddiag_steady(rho(nr_down:nr),number,order,1,0);%下方rho方向的导数系数
d2x_down=d2x0;                          d1x_down=d1x0;
kdx_down=kdx0+nr_down-1;                %极向坐标都是从nt1+1开始，只包含nt_inner个点
%将矩阵转化为正确的格式方便之后在循环中赋值
[d2y0,d1y0,kdy0]=getddiag_steady(theta(nt1+1:nt1+nt_inner),number,order,2,1);
d2y_inner=d2y0;                          d1y_inner=d1y0;
kdy_inner=kdy0+nt1;                      %极向坐标都是从nt1+1开始，只包含nt_inner个点
[d2y0,d1y0,kdy0]=getddiag_steady([theta(t_min:nt1),theta(nt1+nt_inner+1:t_max)-2*pi],number,order,2,0);
d2y_pri=d2y0;                                        d1y_pri=d1y0;
kdy_pri=kdy0+t_min-1+(kdy0>nt1-t_min+1)*nt_inner;    %极向坐标都是从t_min开始，只包含nt_inner个点
[d2y0,d1y0,kdy0]=getddiag_steady(theta(t_min:t_max),number,order,2,0);
d2y_SOL=d2y0;                            d1y_SOL=d1y0;
kdy_SOL=kdy0+t_min-1;                    %极向坐标都是从t_min开始，只包含nt_inner个点

for i=1:nr_inner
    kk=nt1+1:nt1+nt_inner;%本部分表示的是核心区域等离子体
    kdy(i,kk,:) = kdy_inner;
    d2y(i,kk,:) = d2y_inner;
    d1y(i,kk,:) = d1y_inner;
    if i>=nr_down%-- extend t before evaluation 对于私有区域的点
        kk=[t_min:nt1,nt1+nt_inner+1:t_max];
        d2y(i,kk,:)=d2y_pri;     d1y(i,kk,:)=d1y_pri;
        kdy(i,kk,:)=kdy_pri;     
    end
end
for i=nr_inner+1:nr
    kk=(t_min:t_max);
    d2y(i,kk,:)=d2y_SOL;     d1y(i,kk,:)=d1y_SOL;   kdy(i,kk,:)=kdy_SOL;     
end
for j=nt1+1:nt1+nt_inner%本部分表示的是x点上方区域
    d1x(:,j,:) = d1x_up;     d2x(:,j,:) = d2x_up;   kdx(:,j,:) = kdx_up;     
end
for j=[t_min:nt1,nt1+nt_inner+1:t_max]%本部分表示的是x点下方区域
    kk=nr_down:nr;
    d2x(kk,j,:)=d2x_down;    d1x(kk,j,:)=d1x_down;  kdx(kk,j,:)=kdx_down;     
end
%%
xn_range=12;
order_x=order+2;
%在X点附近使用高阶格式而非稳定格式,附近121个点
[d2x_up,d1x_up,kdx_up]=getddiag_steady(rho,number,order_x,1,0);%上方rho方向的导数系数
[d2x0,d1x0,kdx0]=getddiag_steady(rho(nr_down:nr),number,order_x,1,0);%下方rho方向的导数系数
d2x_down=d2x0;                          d1x_down=d1x0;
kdx_down=kdx0+nr_down-1;                %极向坐标都是从nt1+1开始，只包含nt_inner个点
%将矩阵转化为正确的格式方便之后在循环中赋值
[d2y0,d1y0,kdy0]=getddiag_steady(theta(nt1+1:nt1+nt_inner),number,order_x,2,1);
d2y_inner=d2y0;                          d1y_inner=d1y0;
kdy_inner=kdy0+nt1;                      %极向坐标都是从nt1+1开始，只包含nt_inner个点
[d2y0,d1y0,kdy0]=getddiag_steady([theta(t_min:nt1),theta(nt1+nt_inner+1:t_max)-2*pi],number,order_x,2,0);
d2y_pri=d2y0;                                        d1y_pri=d1y0;
kdy_pri=kdy0+t_min-1+(kdy0>nt1-t_min+1)*nt_inner;    %极向坐标都是从t_min开始，只包含nt_inner个点
[d2y0,d1y0,kdy0]=getddiag_steady(theta(t_min:t_max),number,order_x,2,0);
d2y_SOL=d2y0;                            d1y_SOL=d1y0;
kdy_SOL=kdy0+t_min-1;                    %极向坐标都是从t_min开始，只包含nt_inner个点

for i=nr_inner+1-xn_range:nr_inner
    kk=nt1+1:nt1+nt_inner;%本部分表示的是核心区域等离子体
    X_index=[1:xn_range+1,nt_inner+1-xn_range:nt_inner];
    kdy(i,kk(X_index),:) = kdy_inner(1,X_index,:);
    d2y(i,kk(X_index),:) = d2y_inner(1,X_index,:);
    d1y(i,kk(X_index),:) = d1y_inner(1,X_index,:);
    if i>=nr_down%-- extend t before evaluation 对于私有区域的点
        kk=[t_min:nt1,nt1+nt_inner+1:t_max];
        X_index=(nt1-t_min+1+1-xn_range):(nt1-t_min+1+1+xn_range);
        d2y(i,kk(X_index),:)=d2y_pri(1,X_index,:);     d1y(i,kk(X_index),:)=d1y_pri(1,X_index,:);
        kdy(i,kk(X_index),:)=kdy_pri(1,X_index,:);     
    end
end

for i=nr_inner+1:nr_inner+1+xn_range
    kk=(t_min:t_max);
    X_index=[(nt1-t_min+1+1-xn_range):(nt1-t_min+1+1+xn_range) (nt_inner+nt1+1-t_min+1-xn_range):(nt_inner+nt1+1-t_min+1+xn_range)];
    d2y(i,kk(X_index),:)=d2y_SOL(1,X_index,:);     d1y(i,kk(X_index),:)=d1y_SOL(1,X_index,:);   kdy(i,kk(X_index),:)=kdy_SOL(1,X_index,:);     
end
for j=[nt1+1:nt1+1+xn_range,nt1+nt_inner+1-xn_range:nt1+nt_inner]%本部分表示的是x点上方区域
    kk=1:nr;
    X_index=(nr_inner+1-xn_range):(nr_inner+1+xn_range);
    d2x(kk(X_index),j,:) = d2x_up(X_index,1,:);     d1x(kk(X_index),j,:) = d1x_up(X_index,1,:); kdx(kk(X_index),j,:) = kdx_up(X_index,1,:);     
end
for j=  [nt1+1-xn_range:nt1,nt1+nt_inner+1:nt1+nt_inner+1+xn_range]%本部分表示的是x点下方区域
    kk=nr_down:nr;
    X_index=(nr_inner-nr_down+1+1-xn_range):(nr_inner-nr_down+1+1+xn_range);
    d2x(kk(X_index),j,:)=d2x_down(X_index,1,:);     d1x(kk(X_index),j,:)=d1x_down(X_index,1,:); kdx(kk(X_index),j,:)=kdx_down(X_index,1,:);     
end
%%
xn_range=6;
order_x=order+4;
%在X点附近使用高阶格式而非稳定格式,附近121个点
[d2x_up,d1x_up,kdx_up]=getddiag_steady(rho,number,order_x,1,0);%上方rho方向的导数系数
[d2x0,d1x0,kdx0]=getddiag_steady(rho(nr_down:nr),number,order_x,1,0);%下方rho方向的导数系数
d2x_down=d2x0;                          d1x_down=d1x0;
kdx_down=kdx0+nr_down-1;                %极向坐标都是从nt1+1开始，只包含nt_inner个点
%将矩阵转化为正确的格式方便之后在循环中赋值
[d2y0,d1y0,kdy0]=getddiag_steady(theta(nt1+1:nt1+nt_inner),number,order_x,2,1);
d2y_inner=d2y0;                          d1y_inner=d1y0;
kdy_inner=kdy0+nt1;                      %极向坐标都是从nt1+1开始，只包含nt_inner个点
[d2y0,d1y0,kdy0]=getddiag_steady([theta(t_min:nt1),theta(nt1+nt_inner+1:t_max)-2*pi],number,order_x,2,0);
d2y_pri=d2y0;                                        d1y_pri=d1y0;
kdy_pri=kdy0+t_min-1+(kdy0>nt1-t_min+1)*nt_inner;    %极向坐标都是从t_min开始，只包含nt_inner个点
[d2y0,d1y0,kdy0]=getddiag_steady(theta(t_min:t_max),number,order_x,2,0);
d2y_SOL=d2y0;                            d1y_SOL=d1y0;
kdy_SOL=kdy0+t_min-1;                    %极向坐标都是从t_min开始，只包含nt_inner个点

for i=nr_inner+1-xn_range:nr_inner
    kk=nt1+1:nt1+nt_inner;%本部分表示的是核心区域等离子体
    X_index=[1:xn_range+1,nt_inner+1-xn_range:nt_inner];
    kdy(i,kk(X_index),:) = kdy_inner(1,X_index,:);
    d2y(i,kk(X_index),:) = d2y_inner(1,X_index,:);
    d1y(i,kk(X_index),:) = d1y_inner(1,X_index,:);
    if i>=nr_down%-- extend t before evaluation 对于私有区域的点
        kk=[t_min:nt1,nt1+nt_inner+1:t_max];
        X_index=(nt1-t_min+1+1-xn_range):(nt1-t_min+1+1+xn_range);
        d2y(i,kk(X_index),:)=d2y_pri(1,X_index,:);     d1y(i,kk(X_index),:)=d1y_pri(1,X_index,:);
        kdy(i,kk(X_index),:)=kdy_pri(1,X_index,:);     
    end
end

for i=nr_inner+1:nr_inner+1+xn_range
    kk=(t_min:t_max);
    X_index=[(nt1-t_min+1+1-xn_range):(nt1-t_min+1+1+xn_range) (nt_inner+nt1+1-t_min+1-xn_range):(nt_inner+nt1+1-t_min+1+xn_range)];
    d2y(i,kk(X_index),:)=d2y_SOL(1,X_index,:);     d1y(i,kk(X_index),:)=d1y_SOL(1,X_index,:);   kdy(i,kk(X_index),:)=kdy_SOL(1,X_index,:);     
end
for j=[nt1+1:nt1+1+xn_range,nt1+nt_inner+1-xn_range:nt1+nt_inner]%本部分表示的是x点上方区域
    kk=1:nr;
    X_index=(nr_inner+1-xn_range):(nr_inner+1+xn_range);
    d2x(kk(X_index),j,:) = d2x_up(X_index,1,:);     d1x(kk(X_index),j,:) = d1x_up(X_index,1,:); kdx(kk(X_index),j,:) = kdx_up(X_index,1,:);     
end
for j=  [nt1+1-xn_range:nt1,nt1+nt_inner+1:nt1+nt_inner+1+xn_range]%本部分表示的是x点下方区域
    kk=nr_down:nr;
    X_index=(nr_inner-nr_down+1+1-xn_range):(nr_inner-nr_down+1+1+xn_range);
    d2x(kk(X_index),j,:)=d2x_down(X_index,1,:);     d1x(kk(X_index),j,:)=d1x_down(X_index,1,:); kdx(kk(X_index),j,:)=kdx_down(X_index,1,:);     
end
%%
%由于X点处的度规是奇异的，故而该处的度规值不能应用在求导格式中
%下面部分将分界线和十字架的X点处的差分格式更新一下
%在X点附近奇异性明显，应选用高阶格式更好处理
need_change=floor(number/2);
%上方rho方向的导数系数
j=nt1+1;kk=1:nr_inner+1;
[d2x_x,d1x_x,kdx_x]=getddiag_steady(rho(kk),number,order_x,1,0);
to_change=nr_inner+1-need_change:nr_inner+1;
d2x(kk(to_change),j,:)=d2x_x(to_change,:,:);
d1x(kk(to_change),j,:)=d1x_x(to_change,:,:);
kdx(kk(to_change),j,:)=kdx_x(to_change,:,:)+kk(1)-1;
%右侧rho方向的导数系数
j=nt1+1;kk=nr_inner+1:nr;
[d2x_x,d1x_x,kdx_x]=getddiag_steady(rho(kk),number,order_x,1,0);
to_change=1:1+need_change;
d2x(kk(to_change),j,:)=d2x_x(to_change,:,:);
d1x(kk(to_change),j,:)=d1x_x(to_change,:,:);
kdx(kk(to_change),j,:)=kdx_x(to_change,:,:)+kk(1)-1;
%下方rho方向的导数系数
j=nt1+nt_inner+1;kk=nr_down:nr_inner+1;
[d2x_x,d1x_x,kdx_x]=getddiag_steady(rho(kk),number,order_x,1,0);
to_change=nr_inner+1-nr_down+1-need_change:nr_inner+1-nr_down+1;
d2x(kk(to_change),j,:)=d2x_x(to_change,:,:);
d1x(kk(to_change),j,:)=d1x_x(to_change,:,:);
kdx(kk(to_change),j,:)=kdx_x(to_change,:,:)+kk(1)-1;
%左方rho方向的导数系数
j=nt1+nt_inner+1;kk=nr_inner+1:nr;
[d2x_x,d1x_x,kdx_x]=getddiag_steady(rho(kk),number,order_x,1,0);%上方rho方向的导数系数
to_change=1:1+need_change;
d2x(kk(to_change),j,:)=d2x_x(to_change,:,:);
d1x(kk(to_change),j,:)=d1x_x(to_change,:,:);
kdx(kk(to_change),j,:)=kdx_x(to_change,:,:)+kk(1)-1;

%分界面上的特殊极向导数值
i=nr_inner+1;
%右上方分界面
[d2y_x,d1y_x,kdy_x]=getddiag_steady(theta(nt1+1:nt1+nt_inner+1),number,order_x,2,0);
kk=nt1+1:nt1+nt_inner+1;X_index=1:1+need_change;%x点是否更新都可以，因为不需要x点处的差分格式来计算
kdy(i,kk(X_index),:) = kdy_x(1,X_index,:)+kk(1)-1;
d2y(i,kk(X_index),:) = d2y_x(1,X_index,:);
d1y(i,kk(X_index),:) = d1y_x(1,X_index,:);

%左上方分界面
[d2y_x,d1y_x,kdy_x]=getddiag_steady(theta(nt1+1:nt1+nt_inner+1),number,order_x,2,0);
kk=nt1+1:nt1+nt_inner+1;X_index=nt_inner+1-need_change:nt_inner+1;%x点是否更新都可以，因为不需要x点处的差分格式来计算
kdy(i,kk(X_index),:) = kdy_x(1,X_index,:)+kk(1)-1;
d2y(i,kk(X_index),:) = d2y_x(1,X_index,:);
d1y(i,kk(X_index),:) = d1y_x(1,X_index,:);

%左下方分界面
[d2y_x,d1y_x,kdy_x]=getddiag_steady(theta(nt1+nt_inner+1:t_max),number,order_x,2,0);
kk=nt1+nt_inner+1:t_max;X_index=1:1+need_change;%x点是否更新都可以，因为不需要x点处的差分格式来计算
kdy(i,kk(X_index),:) = kdy_x(1,X_index,:)+kk(1)-1;
d2y(i,kk(X_index),:) = d2y_x(1,X_index,:);
d1y(i,kk(X_index),:) = d1y_x(1,X_index,:);

%右下方分界面
[d2y_x,d1y_x,kdy_x]=getddiag_steady(theta(t_min:nt1+1),number,order_x,2,0);
kk=t_min:nt1+1;X_index=(nt1+1-t_min+1)-need_change:(nt1+1-t_min+1);%x点是否更新都可以，因为不需要x点处的差分格式来计算
kdy(i,kk(X_index),:) = kdy_x(1,X_index,:)+kk(1)-1;
d2y(i,kk(X_index),:) = d2y_x(1,X_index,:);
d1y(i,kk(X_index),:) = d1y_x(1,X_index,:);


end



function [d2,d1,kd]=getddiag_steady(x,nd,m,x_or_y,ifperiodic)
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
%     nd, 差分格式中使用的取点数
%     x_or_y:确定此次调用求解的是ρ方向还是θ方向
%     ifperiodic：确定此次调用求解的模块是否为周期性边界条件，=0时表示非周期条件
%*****************************************************
nx=length(x);
ng=floor(nd/2);
kd0=-ng +(0:nd-1);%为确保在均匀格点下使用的是对称格式，必须选择奇数
d2=zeros(nx,nd);
d1=zeros(nx,nd);
kd=zeros(nx,nd);
if ifperiodic==0%非周期条件
    for j=1:nx
        jkd = snearby(j,kd0,nx);%寻找j位置附近的几个点
        d2(j,:)=diff_steady_all(x(j),x(jkd),2,m);
        d1(j,:)=diff_steady_all(x(j),x(jkd),1,m);
        kd(j,:) = jkd;%记录求导数用的哪些点
    end
else
    for j=1:nx
        jkd=mod(kd0+j-1,nx)+1;
        x_choose=x(jkd)+floor((kd0+j-1)/nx)*2*pi;
        d2(j,:)=diff_steady_all(x(j),x_choose,2,m);
        d1(j,:)=diff_steady_all(x(j),x_choose,1,m);
        kd(j,:) = jkd;%记录求导数用的哪些
    end
end
if x_or_y==1%求r方向的差分时
    d2=permute(d2,[1,3,2]);d1=permute(d1,[1,3,2]);kd=permute(kd,[1,3,2]);
else%求t方向的差分时
    d2=permute(d2,[3,1,2]);d1=permute(d1,[3,1,2]);kd=permute(kd,[3,1,2]);
end
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

function c = diff_steady_all(xbar,x,k,m)
%程序要求x以行向量形式输入
%求用x处函数值表达xbar处k阶导数系数,且具有m阶收敛精度，要求x的长度必须大于等于k+m，
%当x的长度为k+m时，退回到普通方法
%输出格式与x相同
L=length(x);
%m=min(m,L-k);%当输入的阶数要求过大时，即退回到普通方法
if L<=k+m
c = fdcoeffF(k,xbar,x);
else
N=k+m;%约束条件的数目
h=x-xbar;
factor=mean(abs(diff(h)));
h=h/factor;%通过近似归一化的操作，使得矩阵的奇异性降低，使得计算更加准确
H=zeros(1,2*N-1);
for i=1:2*N-1
    H(i)=sum(h.^(i-1));%程序中将0^0计算为1
end
A=zeros(N,N);
j=1:N;
for i=1:N
    A(i,:)=H(1,i:i+N-1)./factorial(j-1)/factorial(i-1);
end
B=zeros(N,1);B(k+1)=1;
D=A\B;
M=zeros(N,L);
for i=1:N
    M(i,:)=h.^(i-1)/factorial(i-1);
end
c=D'*M;
c=c/factor^k;
end
end

