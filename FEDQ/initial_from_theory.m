function eq=initial_from_theory(nr,nt)
icheck=0;
%选用较密集的网格保证边界条件足够准确
R1=5;
R2=7;
Z2=4.5;
A=-0.0005;
B=0.01;
rg=linspace(3.8,8,50*25);
zg=linspace(-7,5,150*25);
[R,Z]=meshgrid(rg,zg);
psi=theory_sample(R,Z);
Pprime0=-2./Z2.^2-8./R2.^2;
FFprime0=2*R1^2/Z2^2;
%{
figure
hold on
plot(R(1:25*2:end,1:25*2:end),Z(1:25*2:end,1:25*2:end),'g-')
plot(R(1:25*2:end,1:25*2:end)',Z(1:25*2:end,1:25*2:end)','g-')
contour(R,Z,psi,100);
axis equal
contour(R,Z,psi,[0 0]-0.6,'r.-')
contour(R,Z,psi,[0 0]-1.6,'b.-')
%}
Rlim=[4 10 10 4 4];
%Rlim=[3.9 10 10 3.9 3.9];
Zlim=[-5 -5  2 2 -5];
%Zlim=[-4.9 -4.9  2 2 -4.9];
R0=6.22;
Z0=-1.33;
psi_axis=-4.8;
%[psi_axis,I]=min(psi,[],'all');
%scatter(R(I),Z(I),'o');
%确定一个合适的最外层磁面
psi_lcs=-0.6;
psi_lcs=-0.8;
figure(1)
h=contour(rg,zg,psi,[0 0]+psi_lcs,'Visible','off');
rboundary=h(1,2:1+h(2,1));
zboundary=h(2,2:1+h(2,1));

%plot(rboundary,zboundary)

%确定一个合适的私有区磁面
psi_pri=-1.6;
psi_pri=-1.2;
h=contour(rg,zg,psi,[0 0]+psi_pri,'Visible','off');
rdown=h(1,2:1+h(2,1));
zdown=h(2,2:1+h(2,1));

%plot(rdown,zdown)

close 1
[xup,yup,iiup] = polyxpoly(Rlim,Zlim,rboundary,zboundary);
[xdown,ydown,iidown] = polyxpoly(Rlim,Zlim,rdown,zdown);
divertor_left1=[xup(1) yup(1)];
divertor_left2=[xdown(1) ydown(1)];
divertor_right1=[xup(2) yup(2)];
divertor_right2=[xdown(2) ydown(2)];
divertor_left=divertor_left1+linspace(0,1,100)'*(divertor_left2-divertor_left1);
divertor_right=divertor_right1+linspace(0,1,100)'*(divertor_right2-divertor_right1);%记录顺序均是从上至下
psi_left =interp2(R,Z,psi,divertor_left(:,1),divertor_left(:,2));   psi_left(1)=psi_lcs;    psi_left(end)=psi_pri;
psi_right=interp2(R,Z,psi,divertor_right(:,1),divertor_right(:,2)); psi_right(1)=psi_lcs;   psi_right(end)=psi_pri;
Rdiv_left=divertor_left(:,1)';
Zdiv_left=divertor_left(:,2)';
Rdiv_right=divertor_right(:,1)';
Zdiv_right=divertor_right(:,2)';
Rlcs=[xup(1) rboundary(iiup(1,2)+1:iiup(2,2)) xup(2)];
Zlcs=[yup(1) zboundary(iiup(1,2)+1:iiup(2,2)) yup(2)];
Rpri=[xdown(1) rdown(iidown(1,2)+1:iidown(2,2)) xdown(2)];
Zpri=[ydown(1) zdown(iidown(1,2)+1:iidown(2,2)) ydown(2)];
[Rlcs,Zlcs] = interp2_curve(Rlcs,Zlcs,'linear',1000);
[Rpri,Zpri] = interp2_curve(Rpri,Zpri,'linear',250);
boundary=struct('Rlcs',Rlcs,'Zlcs',Zlcs, ...
    'Rpri',Rpri,'Zpri',Zpri, ...
    'Rdiv_left',Rdiv_left,'Rdiv_right',Rdiv_right,'Zdiv_left',Zdiv_left,'Zdiv_right',Zdiv_right, ...
    'psi_left',psi_left,'psi_right',psi_right,'psi_lcs',psi_lcs,'psi_pri',psi_pri);

%为初始化网格选择较低的网格密度
rg=linspace(4,8,64*10);
zg=linspace(-7,5,128*10);
% if exist("nr","var")
%     rg=linspace(4,8,64*5*floor(nr/64));
%     zg=linspace(-7,5,128*5*floor(nr/64));
% else
%     rg=linspace(4,8,64*5);
%     zg=linspace(-7,5,128*5);
% end

[R,Z]=meshgrid(rg,zg);
psi=theory_sample(R,Z)+0*rand(size(R));
%初始化磁轴位置
icheck=0;
axis_range=0.1;%表示磁轴附近区域的大小，初始化时不在意精度
%[R0,Z0,dpsi_dr,psi_axis,psi_axis_fit]=findaxis(R,Z,rg,zg,psi,psi_lcs,icheck,R0,Z0,psi_axis);
[R0,Z0,dpsi_dr,psi_axis]=fit_axis_initial(R,Z,psi,R0,Z0,axis_range,icheck);
%初始化x点位置和磁通值
icheck=0;
%[Rx,Zx,psi_x,psi] = fit_xpoint(R,Z,psi,psi_lcs,psi_pri,rg,zg,icheck);
Rx_pre=4.6;Zx_pre=-3.9;x_range=0.2;
[Rx,Zx,psi_x,psi] = fit_xpoint_initial(R,Z,psi,psi_lcs,psi_pri,Rx_pre,Zx_pre,x_range,icheck);
%{
figure
hold on
plot(Rdiv_left,Zdiv_left,'k','LineWidth',3);
plot(Rdiv_right,Zdiv_right,'k','LineWidth',3);
plot(Rlcs,Zlcs,'k','LineWidth',3);
plot(Rpri,Zpri,'k','LineWidth',3);
%scatter(Rx,Zx);
%scatter(R0,Z0);
axis equal
axis tight
box off
set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])
axis off
%}
%%

%调节坐标区域范围
Zmax=max(Zlcs);        [~,Zmax_index]=min(abs(zg-Zmax));
Zmin=min(Zpri);        [~,Zmin_index]=min(abs(zg-Zmin));
Rmax=max(Rlcs);        [~,Rmax_index]=min(abs(rg-Rmax));
Rmin=min([Rpri Rlcs]); [~,Rmin_index]=min(abs(rg-Rmin));
index1=max(1,Zmin_index-3):min(Zmax_index+3,length(zg));
index2=max(1,Rmin_index-3):min(Rmax_index+3,length(rg));
psi=psi(index1,index2);
R=R(index1,index2);
Z=Z(index1,index2);
%delta_psi=delta_psi(index1,index2);
rg=rg(index2);
zg=zg(index1);
%%
%将磁通值归一化
psib=psi_lcs-psi_axis;
npsi=(psi-psi_axis)/psib;
ndpsi_dr0=dpsi_dr/psib;
nrho=real(sqrt(npsi));%防止拟合出的磁轴磁通大于最小磁通值所造成的复数问题
npsi_x=(psi_x-psi_axis)/psib;               nrho_x=sqrt(npsi_x);
npsi_pri=(psi_pri-psi_axis)/psib;           nrho_pri=sqrt(npsi_pri);
npsi_left=(psi_left-psi_axis)/psib;         nrho_left=sqrt(npsi_left);  nrho_left(1) =1;    nrho_left(end) =nrho_pri;
npsi_right=(psi_right-psi_axis)/psib;       nrho_right=sqrt(npsi_right);nrho_right(1)=1;    nrho_right(end)=nrho_pri;%用于修正数值计算误差
%%
%旋转坐标系
[alpha,Length]=cart2pol(R0-Rx,Z0-Zx);
[x,y] = rotate(R,Z,R0,Z0,alpha,Length);
[xdiv_left,ydiv_left]=rotate(Rdiv_left,Zdiv_left,R0,Z0,alpha,Length);
[xdiv_right,ydiv_right]=rotate(Rdiv_right,Zdiv_right,R0,Z0,alpha,Length);

%%
%构造极向坐标和径向坐标
if ~exist("nr","var")
nr=64*4;
nt=64*4*4;
end
rho=x_nonuniform_new([0 1],nr,[0 nrho_x],[0.1 0.025],[3 20],0);
%rho=x_nonuniform_new([0 1],nr,[0 nrho_x],[0.1 0.025],[1 10],0);
%rho=x_nonuniform_new([0 1],nr,[0 nrho_x],[0.1 0.025],[3 10],0);
rho = x_adjustment(rho,[0 1],[nrho_pri,nrho_x]);
nr_down =sum(rho<nrho_pri)+1;
nr_inner=sum(rho<nrho_x);
nr_outer=nr-nr_inner;

%theta=x_nonuniform_new([-pi/2 pi/2*5],nt,[0 2*pi 0 2*pi 4.9/4*pi 0 2*pi],[pi/24 pi/24 pi/48/4 pi/48/4 pi/6 pi/9 pi/9],[5 5 40 40 2 2 2]/1.5,0);
%theta=x_nonuniform_new([-pi/2 pi/2*5],nt,[0 2*pi 0 2*pi 4.9/4*pi 0 2*pi],[pi/24 pi/24 pi/48/4 pi/48/4 pi/6 pi/9 pi/9],[5 5 10 10 2 2 2]/1.5,0);
theta=x_nonuniform_new([-pi/2 pi/2*5],nt,[0 2*pi 0 2*pi 4.9/4*pi 0 2*pi],[pi/24 pi/24 pi/48/4 pi/48/4 pi/6 pi/9 pi/9],[5 5 160*2 160*2 2 2 2]/1.5,0);
theta=x_adjustment(theta,[-pi/2 pi/2*5],[0 2*pi]);
nt1=sum(theta<0);
nt2=sum(theta>2*pi)+1;
nt_inner=nt-nt1-nt2;
theta_core=theta(nt1+1:nt-nt2);
theta_pri_right=theta(1:nt1);
theta_pri_left=theta(nt-nt2+1:end);
theta_pri=[theta_pri_right theta_pri_left];

theta_div_left  = coordinate_construct(xdiv_left, ydiv_left);
theta_div_right = coordinate_construct(xdiv_right,ydiv_right);%表示偏滤器上的点对应的θ值，与偏滤器位置记录的点一一对应
theta_div_left_new =interp1(nrho_left, theta_div_left, rho(nr_down:end));
theta_div_right_new=interp1(nrho_right,theta_div_right, rho(nr_down:end));%表示由偏滤器限制的，每个磁面上的θ的边界值
theta_div_left_new =[ones(1,nr_down-1)*pi/2*5  theta_div_left_new];
theta_div_right_new=[ones(1,nr_down-1)*(-pi/2) theta_div_right_new];%将其转化为与所有磁面都有对应值的形式方便下面的计算

t_boundary=zeros(3,nr);
for j=1:nr_down-1
    t_boundary(1,j)=nt1+1;
    t_boundary(2,j)=nt1+nt_inner;
    t_boundary(3,j)=nt_inner;
end
for j=nr_down:nr
    t_boundary(1,j)=sum(theta<theta_div_right_new(j))+1;
    t_boundary(2,j)=sum(theta<theta_div_left_new(j));
    t_boundary(3,j)=t_boundary(2,j)-t_boundary(1,j)+1;
end
t_br1=min(t_boundary(1,nr_down:end));   t_br2=max(t_boundary(1,nr_down:end));
t_bl1=max(t_boundary(2,nr_down:end));   t_bl2=min(t_boundary(2,nr_down:end));
ng_r=3*max((t_br2-t_br1),1);
ng_l=3*max((t_bl1-t_bl2),1);
ng=5;%表示参与重新分配坐标的完整θ方向点数
%ng=min([10,nt1+1-t_br1-5-3,t_bl1-(nt1+nt_inner++1+5)-3]);%表示参与重新分配坐标的完整θ方向点数
nelse=0;%表示追加的额外节点数
tr_min=t_br1-nelse;  tr_max=t_br2+ng_r;   nf_r=tr_max-tr_min; tf_r=(1:nf_r)/nf_r;
tl_max=t_bl1+nelse;  tl_min=t_bl2-ng_l;   nf_l=tl_max-tl_min; tf_l=(1:nf_l)/nf_l;
t_min=tr_min;   t_max=tl_max;%新的极向坐标边界情况
t_boundary(1,1:nr_down-1)=nt1+1;        t_boundary(1,nr_down:end)=t_min;
t_boundary(2,1:nr_down-1)=nt1+nt_inner; t_boundary(2,nr_down:end)=t_max;
t_boundary(3,1:nr_down-1)=nt_inner;     t_boundary(3,nr_down:end)=t_max-t_min+1;
theta_pri_f_all=zeros(nr_inner+1-nr_down+1,t_max-nt_inner-t_min+1);

tf_r=tf_r.^2;
tf_l=tf_l.^2;
for j=nr_down:nr_inner+1
    theta_pri_f=theta_pri;
    theta_pri_f(tr_max-1:-1:tr_min)=tf_r*(theta_div_right_new(j)-theta_pri(tr_min))+theta_pri_f(tr_max-1:-1:tr_min);%这一段点的顺序由靠近y轴向两边
    theta_pri_f(tl_min-nt_inner+1:tl_max-nt_inner)=tf_l*(theta_div_left_new(j)-theta_pri(tl_max-nt_inner))+theta_pri_f(tl_min-nt_inner+1:tl_max-nt_inner);
    theta_pri_f_all(j-nr_down+1,:)=theta_pri_f(t_min:t_max-nt_inner);
end
geometry=struct('nr',nr,'nr_inner',nr_inner,'nr_outer',nr_outer,'nr_down',nr_down, ...
    'nt',nt,'nt_inner',nt_inner,'nt1',nt1,'nt2',nt2,'t_min',t_min,'t_max',t_max,'alpha',alpha,'Length',Length,...
    't_boundary',t_boundary,'tr_max',tr_max,'tl_min',tl_min);%,'r_boundary',r_boundary);

%%
paras.rho_inner=0.02;
paras.n_inner=sum(rho<paras.rho_inner);
%{
figure
hold on
contourf(x,y,npsi,10);
%contourf(x,y,npsi,rho.^2);
%scatter(x,y,'.')
axis equal
%}
%close all
%初始化磁面坐标系网格
x_new=zeros(nr,nt);
y_new=zeros(nr,nt);
R_new=zeros(nr,nt);
Z_new=zeros(nr,nt);
r_new=zeros(nr,nt);
%R_new(1,nt1+1:nt1+nt_inner)=R0;
%Z_new(1,nt1+1:nt1+nt_inner)=Z0;
R_new(1,:)=R0;
Z_new(1,:)=Z0;
figure(1)%contour函数需要一个图窗
for j=2:nr
    %h=contour(x,y,npsi,rho(j)^2*[1 1]);%find surface positions
    h=contour(x,y,npsi,rho(j)^2*[1 1],'visible','off');%find surface positions
    %if j<paras.n_inner
    if j<=6
    r_new(j,nt1+1:nt1+nt_inner)=sqrt(rho(j)^2./(ndpsi_dr0(1)+(ndpsi_dr0(2)-ndpsi_dr0(1))*sin(theta_core+alpha-pi).^2));
    R_new(j,nt1+1:nt1+nt_inner)=r_new(j,nt1+1:nt1+nt_inner).*cos(theta_core+alpha-pi)+R0;
    Z_new(j,nt1+1:nt1+nt_inner)=r_new(j,nt1+1:nt1+nt_inner).*sin(theta_core+alpha-pi)+Z0;
    elseif j<=nr_down-1
        if h(2,1)<length(h)-1%求出的等高线不只有一条
            if h(2,2)<-1%开始位置记录的是下方的
                h(:,1:h(2,1)+1)=[];%将没有用的下方那一段清空
            end
        end
        jx1=h(1,2:h(2,1));%由于是闭合的，因而等高线向量的首尾是相同的
        jy1=h(2,2:h(2,1));
        jtheta1=coordinate_construct(jx1,jy1);
        
        [jtheta1,index]=sort(jtheta1);jx1=jx1(index);jy1=jy1(index);
        curve = csape([ jtheta1 jtheta1(1)+2*pi],[ jx1 jx1(1); jy1 jy1(1)],'period');
        value = fnval(curve,mod(theta_core-jtheta1(1),2*pi)+jtheta1(1));x1=value(1,:);y1=value(2,:);
        
        x_new(j,nt1+1:nt1+nt_inner)=x1;
        y_new(j,nt1+1:nt1+nt_inner)=y1;
        [R_new(j,nt1+1:nt1+nt_inner),Z_new(j,nt1+1:nt1+nt_inner)] = rotate_inverse(x1,y1,R0,Z0,alpha,Length);
    elseif j<=nr_inner
        theta_pri_f=theta_pri_f_all(j-nr_down+1,:);
        if max(h(2,2:h(2,1)+1))>0%（jx1,jy1）与芯部的磁面对应，标号（jx2,jy2）与私有区的磁面对应
            jx1=h(1,2:h(2,1)+1);
            jy1=h(2,2:h(2,1)+1);
            jx2=h(1,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
            jy2=h(2,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
        else
            jx2=h(1,2:h(2,1)+1);
            jy2=h(2,2:h(2,1)+1);
            jx1=h(1,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
            jy1=h(2,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
        end
        if j==nr_down%当j=nr_down时，对应私有区边界线，不需要再次求解
            [jx2,jy2]=rotate(Rpri,Zpri,R0,Z0,alpha,Length);
        end
        jtheta1=coordinate_construct(jx1,jy1);
        jtheta2=coordinate_construct(jx2,jy2);

        [jtheta1,index]=sort(jtheta1);jx1=jx1(index);jy1=jy1(index);
        curve = csape([ jtheta1 jtheta1(1)+2*pi],[ jx1 jx1(1); jy1 jy1(1)],'period');
        value = fnval(curve,mod(theta_core-jtheta1(1),2*pi)+jtheta1(1));x1=value(1,:);y1=value(2,:);
        curve = csape(mod(jtheta2+pi/2,2*pi)-pi/2,[ jx2;jy2]);
        value = fnval(curve,mod(theta_pri_f+pi,2*pi)-pi);x2=value(1,:);y2=value(2,:);

        x_new(j,nt1+1:nt1+nt_inner)=x1;
        y_new(j,nt1+1:nt1+nt_inner)=y1;
        x_new(j,[t_min:nt1,nt1+nt_inner+1:t_max])=x2;
        y_new(j,[t_min:nt1,nt1+nt_inner+1:t_max])=y2;
        [R_new(j,t_min:t_max),Z_new(j,t_min:t_max)] =...
            rotate_inverse(x_new(j,t_min:t_max),y_new(j,t_min:t_max),R0,Z0,alpha,Length);

    elseif j==nr_inner+1
        theta_pri_f=theta_pri_f_all(j-nr_down+1,:);
        if h(2,1)==length(h)-1%求出的等高线只有一条
            jx=h(1,2:h(2,1)+1);
            jy=h(2,2:h(2,1)+1);
            para_pre=0.1:0.1:0.9;
            jx1=jx(jy>-1&abs(jx)>0.025);jx1=[para_pre*jx1(1) jx1 para_pre*jx1(end)];
            jy1=jy(jy>-1&abs(jx)>0.025);jy1=[para_pre*(jy1(1)+1)-1 jy1 para_pre*(jy1(end)+1)-1];
            jx2=jx(jy<-1&abs(jx)>0.05);[jx2,index_x]=sort(jx2);  number_left=sum(jx2<0);
            jy2=jy(jy<-1&abs(jx)>0.05);jy2=jy2(index_x);
            jx2=[jx2(1:number_left) para_pre*jx2(number_left)       0   para_pre*jx2(number_left+1)         jx2(number_left+1:end)];
            jy2=[jy2(1:number_left) para_pre*(jy2(number_left)+1)-1 -1  para_pre*(jy2(number_left+1)+1)-1   jy2(number_left+1:end)];

            jtheta1=coordinate_construct(jx1,jy1);
            [jtheta1,index_theta]=sort(jtheta1);jx1=jx1(index_theta);jy1=jy1(index_theta);
            jr1=sqrt(jx1.^2+jy1.^2);
            r1=interp1([0 jtheta1 2*pi],[1 jr1 1],theta_core);
            [x1,y1] =  xy_theta_r_new(theta_core,r1);

            jtheta2=coordinate_construct(jx2,jy2);jtheta2(number_left+length(para_pre)+1)=0;
            jr2=sqrt(jx2.^2+(jy2+2).^2);
            r2=interp1(mod(jtheta2+pi/2,2*pi)-pi/2,jr2,mod(theta_pri_f+pi/2,2*pi)-pi/2);
            [x2,y2] =  xy_theta_r_new(theta_pri_f,r2);

            x_new(j,nt1+1:nt1+nt_inner)=x1;
            y_new(j,nt1+1:nt1+nt_inner)=y1;
            x_new(j,[t_min:nt1,nt1+nt_inner+1:t_max])=x2;
            y_new(j,[t_min:nt1,nt1+nt_inner+1:t_max])=y2;
            [R_new(j,t_min:t_max),Z_new(j,t_min:t_max)] =...
                rotate_inverse(x_new(j,t_min:t_max),y_new(j,t_min:t_max),R0,Z0,alpha,Length);
        else
            %求出的等高线有两条时的处理方法
            if max(h(2,2:h(2,1)+1))>0%标号1与芯部的磁面对应，标号2与私有区的磁面对应
                jx_up=h(1,2:h(2,1)+1);
                jy_up=h(2,2:h(2,1)+1);
                jx_down=h(1,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
                jy_down=h(2,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
            else
                jx_down=h(1,2:h(2,1)+1);
                jy_down=h(2,2:h(2,1)+1);
                jx_up=h(1,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
                jy_up=h(2,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
            end
            para_pre=0.1:0.1:0.9;
            jx1=jx_up(jy_up>-1&abs(jx_up)>0.025);jx1=[ para_pre*jx1(1) jx1 para_pre*jx1(end) ];
            jy1=jy_up(jy_up>-1&abs(jx_up)>0.025);jy1=[ para_pre*(jy1(1)+1)-1 jy1 para_pre*(jy1(end)+1)-1 ];
            jx2=jx_down(jy_down<-1&abs(jx_down)>0.05);[jx2,index_x]=sort(jx2);  number_left=sum(jx2<0);
            jy2=jy_down(jy_down<-1&abs(jx_down)>0.05);jy2=jy2(index_x);
            jx2=[jx2(1:number_left) para_pre*jx2(number_left)       0   para_pre*jx2(number_left+1)         jx2(number_left+1:end)];
            jy2=[jy2(1:number_left) para_pre*(jy2(number_left)+1)-1 -1  para_pre*(jy2(number_left+1)+1)-1   jy2(number_left+1:end)];

            jtheta1=coordinate_construct(jx1,jy1);[jtheta1,index_theta]=sort(jtheta1);jx1=jx1(index_theta);jy1=jy1(index_theta);
            jtheta2=coordinate_construct(jx2,jy2);jtheta2(number_left+length(para_pre)+1)=0;
            index_right=number_left+length(para_pre)+1:length(jtheta2);index_left=1:number_left+length(para_pre)+1;
            
            curve1 = csape([0 jtheta1 2*pi],[0 jx1 0;-1 jy1 -1]);
            value1 = fnval(curve1,mod(theta_core-jtheta1(1),2*pi)+jtheta1(1));x1=value1(1,:);y1=value1(2,:);
            index_right_new=1:nt1-t_min+1;index_left_new=nt1-t_min+2:t_max-t_min+1-nt_inner;
            curve2_right = csape(mod(jtheta2(index_right)+pi,2*pi)-pi,[jx2(index_right);jy2(index_right)]);
            curve2_left =  csape(mod(jtheta2(index_left) +pi,2*pi)-pi,[jx2(index_left );jy2(index_left )]);
            value2_right= fnval(curve2_right,mod(theta_pri_f(index_right_new)+pi,2*pi)-pi);x2_right=value2_right(1,:);y2_right=value2_right(2,:);
            value2_left = fnval(curve2_left ,mod(theta_pri_f(index_left_new )+pi,2*pi)-pi);x2_left= value2_left(1,:); y2_left =value2_left (2,:);
            x2=[x2_right x2_left];y2=[y2_right y2_left];

            x_new(j,nt1+1:nt1+nt_inner)=x1;                 x_new(j,nt1+1)=0;
            y_new(j,nt1+1:nt1+nt_inner)=y1;                 y_new(j,nt1+1)=-1;
            x_new(j,[t_min:nt1,nt1+nt_inner+1:t_max])=x2;   x_new(j,nt1+nt_inner+1)=0;
            y_new(j,[t_min:nt1,nt1+nt_inner+1:t_max])=y2;   y_new(j,nt1+nt_inner+1)=-1;

            [R_new(j,t_min:t_max),Z_new(j,t_min:t_max)] =...
                rotate_inverse(x_new(j,t_min:t_max),y_new(j,t_min:t_max),R0,Z0,alpha,Length);
        end
    elseif j<=nr%对于SOL区域不能采用θ，r坐标来插值，太不准确，应使用θ，x坐标来插值
        theta_f=theta;
        theta_f(tr_max-1:-1:tr_min)=tf_r*(theta_div_right_new(j)-theta(tr_min))+theta_f(tr_max-1:-1:tr_min);%这一段点的顺序由靠近y轴向两边
        theta_f(tl_min+1:tl_max)=tf_l*(theta_div_left_new(j)-theta(tl_max))+theta_f(tl_min+1:tl_max);
        theta_f=theta_f(t_min:t_max);
        if j==nr%此处为处理最外层边界位置
            [jx1,jy1]=rotate(Rlcs,Zlcs,R0,Z0,alpha,Length);
        else
            jx1=h(1,2:h(2,1)+1);
            jy1=h(2,2:h(2,1)+1);
        end
        jtheta1=coordinate_construct(jx1,jy1);
        curve = csape(jtheta1,[ jx1;jy1]);
        %{
        figure(11)
        hold on
        plot(jx1,jy1,'r*-');
        fnplt(curve,'b.-');
        %}
        value = fnval(curve,theta_f);x1=value(1,:);y1=value(2,:);
        x_new(j,t_min:t_max)=x1;y_new(j,t_min:t_max)=y1;
        [R_new(j,t_min:t_max),Z_new(j,t_min:t_max)] =...
            rotate_inverse(x_new(j,t_min:t_max),y_new(j,t_min:t_max),R0,Z0,alpha,Length);
    else
    end
    if j<nr_down
        R_new(j,1:nt1+1-1)=R_new(j,nt1+1);  Z_new(j,1:nt1+1-1)=Z_new(j,nt1+1);
        R_new(j,nt1+nt_inner+1:end)=R_new(j,nt1+nt_inner);Z_new(j,nt1+nt_inner+1:end)=Z_new(j,nt1+nt_inner);
    else
        R_new(j,1:t_min-1)=R_new(j,t_min);  Z_new(j,1:t_min-1)=Z_new(j,t_min);
        R_new(j,t_max+1:end)=R_new(j,t_max);Z_new(j,t_max+1:end)=Z_new(j,t_max);
    end
end
close 1
[Theta,Rho]=meshgrid(theta,rho);
npsi=Rho.^2;
psi=npsi*psib+psi_axis;
%%
%{
index=false(nr,nt);
index(5:nr_inner,[nt1-5:nt1+5 nt1+nt_inner-5:nt1+nt_inner+5])=true;
x_new(index)=x_theta_y(Theta(index),y_new(index));
[R_new(index),Z_new(index)]=rotate_inverse(x_new(index),y_new(index),R0,Z0,alpha,Length);
%将芯部和私有区极向坐标最密集的区域保证按θ划分
index=false(nr,nt);
index(nr_inner+2:nr,[nt1-5:nt1+5 nt1+nt_inner-5:nt1+nt_inner+5])=true;
y_new(index)=y_theta_x(Theta(index),x_new(index));
[R_new(index),Z_new(index)]=rotate_inverse(x_new(index),y_new(index),R0,Z0,alpha,Length);


r(1:nr,nt1+1:nt1+nt_inner)=sqrt(x_new(1:nr,nt1+1:nt1+nt_inner).^2+(y_new(1:nr,nt1+1:nt1+nt_inner)).^2);
r(nr_down:nr,[ng+5+t_min:nt1 nt1+nt_inner+1:t_max-ng-5])=...
    sqrt(x_new(nr_down:nr,[ng+5+t_min:nt1 nt1+nt_inner+1:t_max-ng-5]).^2+(y_new(nr_down:nr,[ng+5+t_min:nt1 nt1+nt_inner+1:t_max-ng-5])+2).^2);
[x_new(1:nr,nt1+1:nt1+nt_inner),y_new(1:nr,nt1+1:nt1+nt_inner)] = ...
    xy_theta_r_new(Theta(1:nr,nt1+1:nt1+nt_inner),r(1:nr,nt1+1:nt1+nt_inner),1);
[x_new(nr_down:nr,[ng+5+t_min:nt1 nt1+nt_inner+1:t_max-ng-5]),y_new(nr_down:nr,[ng+5+t_min:nt1 nt1+nt_inner+1:t_max-ng-5])] = ...
    xy_theta_r_new(Theta(nr_down:nr,[ng+5+t_min:nt1 nt1+nt_inner+1:t_max-ng-5]),r(nr_down:nr,[ng+5+t_min:nt1 nt1+nt_inner+1:t_max-ng-5]),-1);
[R_new(1:nr,nt1+1:nt1+nt_inner),Z_new(1:nr,nt1+1:nt1+nt_inner)]=...
    rotate_inverse(x_new(1:nr,nt1+1:nt1+nt_inner),y_new(1:nr,nt1+1:nt1+nt_inner),R0,Z0,alpha,Length);
[R_new(nr_down:nr,[ng+5+t_min:nt1 nt1+nt_inner+1:t_max-ng-5]),Z_new(nr_down:nr,[ng+5+t_min:nt1 nt1+nt_inner+1:t_max-ng-5])]=...
    rotate_inverse(x_new(nr_down:nr,[ng+5+t_min:nt1 nt1+nt_inner+1:t_max-ng-5]),y_new(nr_down:nr,[ng+5+t_min:nt1 nt1+nt_inner+1:t_max-ng-5]),R0,Z0,alpha,Length);

%}
%%
%{
%本部分用于检测初始化之后的误差
psi_theory = theory_sample(R_new,Z_new);
delta=abs(psi-psi_theory)/psib+realmin;
plot_matrix(R_new,Z_new,log10(delta),'delta',geometry,boundary,[0:-0.5:-6]);
colorbar
axis equal

psi_theory = theory_sample(R_new,Z_new);
delta=abs(psi-psi_theory)/psib+realmin;
plot_matrix(R_new,Z_new,log10(delta),'delta',geometry,boundary,[0:-0.5:-6]);
colorbar
axis equal
%}
%plot_contour(R_new,Z_new,Rx,Zx,geometry,boundary)

%%
Pprime=zeros(size(rho))+Pprime0;
FFprime=zeros(size(rho))+FFprime0;
Pprime_pri=zeros(size(rho))+Pprime0;
FFprime_pri=zeros(size(rho))+FFprime0;
gpsi=rho.^2*psib+psi_axis;
gnpsi=rho.^2;
gpsi_pri=gpsi;
%%
%设置其他参数
paras.imax=40;
paras.err=1e-15;
%所有控制X点范围的参数必须与网格密度挂钩，否则无法在网格加密时提高X点处计算精度
paras.X_y_range=(y_new(nr_inner+1-3,nt1+1)+1)*Length;
paras.X_x_range=(x_new(nr_inner+1+3,nt1+1))*Length;%用于在更新网格时确定需要借助拟合更新的范围
paras.xn_range=2;%分界面两侧各xn_range个面使用拟合方法修改
X_n_fit=5;
paras.theta_range=(theta(nt1+1+3)-theta(nt1+1));%用于在更新网格时确定拟合磁面的权重，使得X点处磁面更准确
paras.rho_range=rho(nr_inner+1)-rho(nr_inner+1-3);
%{
%该部分代码用于设置合适的x点范围控制相关参数，显示x点限制的区域适宜
rx=linspace(0,1.5,100);
tx=linspace(0,2*pi,100);
[x1,y1]=xy_theta_r_new(ones(1,100)*paras.theta_range,rx,1);
[x2,y2]=xy_theta_r_new(2*pi-ones(1,100)*paras.theta_range,rx,1);
[x3,y3]=xy_theta_r_new(-ones(1,100)*paras.theta_range,rx,-1);
[x4,y4]=xy_theta_r_new(2*pi+ones(1,100)*paras.theta_range,rx,-1);
plot_contour(x_new,y_new,0,-1,geometry)
hold on
plot(x1,y1,x2,y2,x3,y3,x4,y4)
%}

eq = struct('R',R_new,'Z',Z_new,'theta',theta,'rho',rho,'geometry',geometry,'boundary',boundary, ...
    'R0',R0,'Z0',Z0,'Rx',Rx,'Zx',Zx,'paras',paras, ...
    'psi',psi,'psi_lcs',psi_lcs,'psib',psib,'psi_axis',psi_axis,'psi_pri',psi_pri,'psi_x',psi_x,...
    'ndpsi_dr',ndpsi_dr0,'npsi',npsi,'npsi_x',npsi_x,'npsi_left',npsi_left,'npsi_right',npsi_right,'npsi_pri',npsi_pri, ...
    'gpsi',gpsi,'gnpsi',gnpsi,'Pprime',Pprime,'FFprime',FFprime,'gpsi_pri',gpsi_pri,'Pprime_pri',Pprime_pri,'FFprime_pri',FFprime_pri);
% eq.M=Metric_SOL_steady(R_new,Z_new,rho,theta,geometry,7,4,X_n_fit,icheck);
% eq.Js=evalJs(gpsi,Pprime,FFprime,gpsi,Pprime_pri,FFprime_pri,R_new,psi,geometry);


end









function [rg,zg,psi,Rlim,Zlim,R0,Z0,psi_axis,eq_right]=data
load('D:\董天佐\毕业设计\MPU_SOL\data.mat');
psi=paras.eqdata.psirz;
rg=paras.eqdata.rg;
zg=paras.eqdata.zg;
R0=paras.R0;
Z0=paras.Z0;
psi_axis=paras.eqdata.psiaxis;
Rlim=paras.eqdata.Rlim;
Zlim=paras.eqdata.Zlim;
Pprime=paras.Pprime;
FFprime=paras.FFprime;
gpsi=paras.gpsi;
psiedge=paras.eqdata.psiedge;
gnpsi=gpsi;
gpsi=gpsi*(psiedge-psi_axis)+psi_axis;
eq_right=struct('Pprime',Pprime,'FFprime',FFprime,'gpsi',gpsi,'gnpsi',gnpsi);
end

function geometry=geometry_input
%geometry=geometry_input;
geometry.nr_inner=2^7;
geometry.nr_outer=16;
geometry.nr_down=16;
geometry.nr=geometry.nr_inner+geometry.nr_outer;
geometry.nt_inner=2^6;
geometry.nt1=geometry.nr_outer+geometry.nr_down+16;
geometry.nt2=geometry.nr_outer+geometry.nr_down+16;
geometry.nt=geometry.nt_inner+geometry.nt1+geometry.nt2;
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


