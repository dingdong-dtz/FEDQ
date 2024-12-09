function [R_new,Z_new]=smooth_1Dxin(R,Z,rho,theta,geometry,rho_x,icheck)

%geometry=eq.geometry;
%rho=eq.rho;
%theta=eq.theta;
%[Rho,Theta]=meshgrid(eq.rho,eq.theta);
%R=eq.R;
%Z=eq.Z;
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

[Theta,Rho]=meshgrid(theta,rho);
Z_new=Z;
R_new=R;
weight=exp(-rho.^2/0.2^2);
weight(1)=100;
index1=rho<2;%内部的网格，估计为线性关系
for j=nt1+1:nt1+nt_inner
    sfR=fit(Rho(index1,j),R(index1,j),'poly2','weight',weight(index1));
    R_new(index1,j)=(sfR(Rho(index1,j))-R(index1,j)).*exp(-Rho(index1,j).^2/0.2^2)+R(index1,j);
    sfZ=fit(Rho(index1,j),Z(index1,j),'poly2','weight',weight(index1));
    Z_new(index1,j)=(sfZ(Rho(index1,j))-Z(index1,j)).*exp(-Rho(index1,j).^2/0.2^2)+Z(index1,j);
    R(index1,j)=R_new(index1,j);
    Z(index1,j)=Z_new(index1,j);
end
%{
weight=ones(size(rho));
weight(1)=100;
weight(nr_inner+1)=100;
weight(end)=100;
for j=nt1+1:nt1+nt_inner
    index1=nr_inner+1:nr;%内部的网格，估计为线性关系
    sfR=fit(Rho(index1,j)-Rho(nr_inner+1,j),R(index1,j),'poly4','weight',weight(index1));
    sfZ=fit(Rho(index1,j)-Rho(nr_inner+1,j),Z(index1,j),'poly4','weight',weight(index1));
    R_new(index1,j)=sfR(Rho(index1,j)-Rho(nr_inner+1,j));
    Z_new(index1,j)=sfZ(Rho(index1,j)-Rho(nr_inner+1,j));
    R(index1,j)=R_new(index1,j);
    Z(index1,j)=Z_new(index1,j);
end
%}
%%
%下面对网格做径向平滑
ww=10000000;
pp=0.999995;
theta_range=pi/144;
rho_range=0.01;
for j=t_min:t_max
    if j>nt1+1&&j<=nt1+nt_inner
        weight=ones(size(rho));
        weight(1)=ww;
        weight(end)=ww;
        weight=weight+ww*exp(-(rho-rho_x).^2/rho_range^2)*exp(-(mod(theta(j)+pi,2*pi)-pi)^2/theta_range^2);
        
        p=ones(size(rho))*pp;
        p=p+(1-pp)*exp(-(rho-rho_x).^2/rho_range^2)*exp(-(mod(theta(j)+pi,2*pi)-pi)^2/theta_range^2);
        value=csaps(rho,[R(:,j) Z(:,j)]',p,rho,weight);
        R_new(:,j)=value(1,:)'; Z_new(:,j)=value(2,:)';
    elseif j==nt1+1
        weight=ones(size(rho));
        weight(1)=ww;
        weight(nr_down)=ww;
        weight(end)=ww;
        weight=weight+100000*exp(-(rho-rho_x).^2/rho_range^2)*exp(-(mod(theta(j)+pi,2*pi)-pi)^2/theta_range^2);
        
        p=ones(size(rho))*pp;
        p=p+(1-pp)*exp(-(rho-rho_x).^2/rho_range^2)*exp(-(mod(theta(j)+pi,2*pi)-pi)^2/theta_range^2);
        value=csaps(rho(1:nr_inner+1),[R(1:nr_inner+1,j) Z(1:nr_inner+1,j)]',p(1:nr_inner+1),rho(1:nr_inner+1),weight(1:nr_inner+1));
        R_new(1:nr_inner+1,j)=value(1,:)'; Z_new(1:nr_inner+1,j)=value(2,:)';
        value=csaps(rho(nr_inner+1:end),[R(nr_inner+1:end,j) Z(nr_inner+1:end,j)]',p(nr_inner+1:end),rho(nr_inner+1:end),weight(nr_inner+1:end));
        R_new(nr_inner+1:end,j)=value(1,:)'; Z_new(nr_inner+1:end,j)=value(2,:)';
    elseif j<=nt1||j>nt1+nt_inner+1
        
        weight=ones(size(rho));
        weight(1)=ww;
        weight(nr_down)=ww;
        weight(end)=ww;
        weight=weight+100000*exp(-(rho-rho_x).^2/rho_range^2)*exp(-(mod(theta(j)+pi,2*pi)-pi)^2/theta_range^2);
        
        p=ones(size(rho))*pp;
        p=p+(1-pp)*exp(-(rho-rho_x).^2/rho_range^2)*exp(-(mod(theta(j)+pi,2*pi)-pi)^2/theta_range^2);
        value=csaps(rho(nr_down:nr),[R(nr_down:nr,j) Z(nr_down:nr,j)]',p(nr_down:nr),rho(nr_down:nr),weight(nr_down:nr));
        R_new(nr_down:nr,j)=value(1,:)'; Z_new(nr_down:nr,j)=value(2,:)';
    elseif j==nt1+nt_inner+1
        weight=ones(size(rho));
        weight(1)=ww;
        weight(nr_down)=ww;
        weight(end)=ww;
        weight=weight+100000*exp(-(rho-rho_x).^2/rho_range^2)*exp(-(mod(theta(j)+pi,2*pi)-pi)^2/theta_range^2);
        
        p=ones(size(rho))*pp;
        p=p+(1-pp)*exp(-(rho-rho_x).^2/rho_range^2)*exp(-(mod(theta(j)+pi,2*pi)-pi)^2/theta_range^2);
        value=csaps(rho(nr_down:nr_inner+1),[R(nr_down:nr_inner+1,j) Z(nr_down:nr_inner+1,j)]',p(nr_down:nr_inner+1),rho(nr_down:nr_inner+1),weight(nr_down:nr_inner+1));
        R_new(nr_down:nr_inner+1,j)=value(1,:)'; Z_new(nr_down:nr_inner+1,j)=value(2,:)';
        value=csaps(rho(nr_inner+1:end),[R(nr_inner+1:end,j) Z(nr_inner+1:end,j)]',p(nr_inner+1:end),rho(nr_inner+1:end),weight(nr_inner+1:end));
        R_new(nr_inner+1:end,j)=value(1,:)'; Z_new(nr_inner+1:end,j)=value(2,:)';
    end
end
R_new(1,:)=R(1,:);Z_new(1,:)=Z(1,:);%磁轴处的位置不应随平滑而变化

if icheck
figure
hold on
f=plot(R(:,nt1+1:nt1+nt_inner),Z(:,nt1+1:nt1+nt_inner),'b.-');
plot(R(nr_down:nr,[t_min:nt1,nt1+nt_inner+1:t_max]),Z(nr_down:nr,[t_min:nt1,nt1+nt_inner+1:t_max]),'b.-')
plot(R(1:nr_inner,[nt1+1:nt1+nt_inner,nt1+1])',Z(1:nr_inner,[nt1+1:nt1+nt_inner,nt1+1])','b.-')
plot(R(nr_down:nr_inner,[t_min:nt1,nt1+nt_inner+1:t_max])',Z(nr_down:nr_inner,[t_min:nt1,nt1+nt_inner+1:t_max])','b.-')
plot(R(1+nr_inner:nr,t_min:t_max)',Z(1+nr_inner:nr,t_min:t_max)','b.-')
%figure
hold on
g=plot(R_new(:,nt1+1:nt1+nt_inner),Z_new(:,nt1+1:nt1+nt_inner),'r.-');
plot(R_new(nr_down:nr,[t_min:nt1,nt1+nt_inner+1:t_max]),Z_new(nr_down:nr,[t_min:nt1,nt1+nt_inner+1:t_max]),'r.-')
plot(R_new(1:nr_inner,[nt1+1:nt1+nt_inner,nt1+1])',Z_new(1:nr_inner,[nt1+1:nt1+nt_inner,nt1+1])','r.-')
plot(R_new(nr_down:nr_inner,[t_min:nt1,nt1+nt_inner+1:t_max])',Z_new(nr_down:nr_inner,[t_min:nt1,nt1+nt_inner+1:t_max])','r.-')
plot(R_new(1+nr_inner:nr,t_min:t_max)',Z_new(1+nr_inner:nr,t_min:t_max)','r.-')
legend([f(1,1);g(1,1)],'before smooth','after smooth')
end
%%
%{
%下面对网格做极向平滑
Z_new_new=Z_new;
R_new_new=R_new;
for i=6:nr
    if i<nr_inner
        index=nt1+1:nt1+nt_inner;
        theta_core=[theta(index)-2*pi theta(index) theta(index)+2*pi];
        R_temp=smooth(theta_core,R_new(i,[index index index]),7,'sgolay',5);
        Z_temp=smooth(theta_core,Z_new(i,[index index index]),7,'sgolay',5);
        R_new_new(i,nt1+1:nt1+nt_inner)=R_temp(nt_inner+1:nt_inner*2);
        Z_new_new(i,nt1+1:nt1+nt_inner)=Z_temp(nt_inner+1:nt_inner*2);
    end
    if i<nr_inner&&i>=nr_down
        index=[t_min:nt1 nt1+nt_inner+1:t_max];
        R_new_new(i,index)=smooth(theta(index),R_new(i,index),7,'sgolay',5);
        Z_new_new(i,index)=smooth(theta(index),Z_new(i,index),7,'sgolay',5);
    end
    if i>=nr_inner+3
        index=[t_min:t_max];
        R_new_new(i,index)=smooth(theta(index),R_new(i,index),7,'sgolay',5);
        Z_new_new(i,index)=smooth(theta(index),Z_new(i,index),7,'sgolay',5);
   end
end
R_new=R_new_new;Z_new=Z_new_new;
if icheck
figure
hold on
f=plot(R(:,nt1+1:nt1+nt_inner),Z(:,nt1+1:nt1+nt_inner),'b.-');
plot(R(nr_down:nr,[t_min:nt1,nt1+nt_inner+1:t_max]),Z(nr_down:nr,[t_min:nt1,nt1+nt_inner+1:t_max]),'b.-')
plot(R(1:nr_inner,[nt1+1:nt1+nt_inner,nt1+1])',Z(1:nr_inner,[nt1+1:nt1+nt_inner,nt1+1])','b.-')
plot(R(nr_down:nr_inner,[t_min:nt1,nt1+nt_inner+1:t_max])',Z(nr_down:nr_inner,[t_min:nt1,nt1+nt_inner+1:t_max])','b.-')
plot(R(1+nr_inner:nr,t_min:t_max)',Z(1+nr_inner:nr,t_min:t_max)','b.-')

g=plot(R_new(:,nt1+1:nt1+nt_inner),Z_new(:,nt1+1:nt1+nt_inner),'r.-');
plot(R_new(nr_down:nr,[t_min:nt1,nt1+nt_inner+1:t_max]),Z_new(nr_down:nr,[t_min:nt1,nt1+nt_inner+1:t_max]),'r.-')
plot(R_new(1:nr_inner,[nt1+1:nt1+nt_inner,nt1+1])',Z_new(1:nr_inner,[nt1+1:nt1+nt_inner,nt1+1])','r.-')
plot(R_new(nr_down:nr_inner,[t_min:nt1,nt1+nt_inner+1:t_max])',Z_new(nr_down:nr_inner,[t_min:nt1,nt1+nt_inner+1:t_max])','r.-')
plot(R_new(1+nr_inner:nr,t_min:t_max)',Z_new(1+nr_inner:nr,t_min:t_max)','r.-')
legend([f(1,1);g(1,1)],'before smooth','after smooth')
end
%}
end

%{
%index=find(theta>=pi/2&theta<=pi/2*3);
index=nt1+1+1:nt1+nt_inner;
for j=index
    weight=ones(size(rho));
    weight(1)=10000000;
    weight(end)=10000000;
    weight=weight+100000*exp(-(rho-rho_x).^2/rho_range^2)*exp(-(mod(theta(j)+pi,2*pi)-pi)^2/theta_range^2);
    pp=0.99999;
    p=ones(size(rho))*pp;
    p=p+(1-pp)*exp(-(rho-rho_x).^2/rho_range^2)*exp(-(mod(theta(j)+pi,2*pi)-pi)^2/theta_range^2);
    %R_new(:,j)=csaps(rho',R(:,j),0.9999,rho',weight);
    %Z_new(:,j)=csaps(rho',Z(:,j),0.9999,rho',weight);
    %value=csaps(rho,[R(:,j) Z(:,j)]',0.9999,rho,weight);
    value=csaps(rho,[R(:,j) Z(:,j)]',p,rho,weight);
    R_new(:,j)=value(1,:)'; Z_new(:,j)=value(2,:)';
end
figure
hold on
plot(R(:,index),Z(:,index),'g.-')
plot(R_new(:,index),Z_new(:,index),'b.-')
plot(R(:,index)',Z(:,index)','g.-')
plot(R_new(:,index)',Z_new(:,index)','b.-')
plot(R(nr_inner+1,index)',Z(nr_inner+1,index)','g.-','LineWidth',2)
plot(R_new(nr_inner+1,index)',Z_new(nr_inner+1,index)','b.-','LineWidth',2)
%}
%{
%在二维空间中直接平滑网格无法很好的近似x点附近的情况
index1=1:nr_inner+1-20;
index2=nt1+1:nt1+nt_inner;
[Theta,Rho]=meshgrid(theta,rho);
Rho_core=Rho(index1,index2);      Rho_core2=[Rho_core Rho_core Rho_core];
Theta_core=Theta(index1,index2);  Theta_core2=[Theta_core-2*pi Theta_core Theta_core+2*pi];
R_core=R(index1,index2);      R_core2=[R_core R_core R_core];
Z_core=Z(index1,index2);      Z_core2=[Z_core Z_core Z_core];
y=[permute(R_core2,[3 1 2]);permute(Z_core2,[3 1 2])];
x = {rho(index1),[theta(index2)-2*pi theta(index2) theta(index2)+2*pi]};
weight=ones(size(rho(index1)));
weight(1)=100000;
weight(end)=100000;
Weight={weight,ones(1,length(index2)*3)};
pp=0.99999;
p=ones(size(rho(index1)))*pp;
p(end-5:end)=1;
P={p,ones(1,length(index2)*3)*pp};
xx = {rho(index1),theta(index2)};
%[value,p]=csaps(x,y,0.9999,xx,Weight);
[value,p]=csaps(x,y,P,xx,Weight);
figure
hold on
plot(R(index1,index2),Z(index1,index2),'g.-')
plot(R(index1,index2)',Z(index1,index2)','g.-')
plot(shiftdim(value(1,:,:)),shiftdim(value(2,:,:)),'r.-')
plot(shiftdim(value(1,:,:))',shiftdim(value(2,:,:))','r.-')
%}