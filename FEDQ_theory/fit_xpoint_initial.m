function [Rx,Zx,psi_x,psi_fit] = fit_xpoint_initial(R,Z,psi,psi_lcs,psi_pri,Rx_pre,Zx_pre,x_range,icheck)
%[Rx,Zx,psi_x,psi_fit] = fit_xpoint_initial(R,Z,psi,psi_lcs,psi_pri,Rx_pre,Zx_pre,x_range,icheck)
%对x点做拟合，找到x点位置和磁通值
%该函数与fit_xpoint的区别为在矩形网格上进行初始化
if ~exist("Rx_pre","var")
[~,Rx_pre,Zx_pre]=estimate_x(psi,psi_lcs,psi_pri,R,Z);
end
index=((R-Rx_pre).^2+(Z-Zx_pre).^2)<x_range^2;
while sum(index(:))<25
x_range=x_range*1.5;
index=((R-Rx_pre).^2+(Z-Zx_pre).^2)<x_range^2;
end
near_R=R(index);
near_Z=Z(index);
near_psi=psi(index);

sf = fit([near_R(:), near_Z(:)],near_psi(:),'poly22');
M=[sf.p20,sf.p11/2;
    sf.p11/2,sf.p02];
T=[sf.p10 sf.p01];
[V,D]=eig(M);%M*V = V*D
e1=V(:,1)/(V(:,1)'*V(:,1));
e2=V(:,2)/(V(:,2)'*V(:,2));
H=[e1,e2];
TH=T*H;
psi_x=sf.p00-TH(1)^2/4/D(1,1)-TH(2)^2/4/D(2,2);
X_Point=H*[-TH(1)/2/D(1,1);-TH(2)/2/D(2,2)];
Rx=X_Point(1);
Zx=X_Point(2);
psi_fit_diff=(sf(R,Z)-psi).*exp(-((R-Rx).^2+(Z-Zx).^2)/x_range^2);
psi_fit=psi+psi_fit_diff;

if icheck
    figure('Name','psi_fit','NumberTitle','off')
    hold on
    contour(R,Z,psi.*index);
    contour(R,Z,psi_fit.*index);
    scatter(Rx_pre,Zx_pre);
    scatter(Rx,Zx);
    axis([Rx_pre-x_range,Rx_pre+x_range,Zx_pre-x_range,Zx_pre+x_range]);
    %{
    near_psi_fit=sf(near_R,near_Z);
    figure('Name','psi_fit','NumberTitle','off')
    hold on
    contour(near_R,near_Z,near_psi);
    contour(near_R,near_Z,near_psi_fit);
    scatter(Rx_pre,Zx_pre);
    scatter(Rx,Zx);
    figure
    hold on
    contour(R,Z,psi,linspace(psi_pri,psi_lcs,25));
    contour(R,Z,psi_fit,linspace(psi_pri,psi_lcs,25));
    contour(R,Z,psi,psi_x*[1 1]);
    contour(R,Z,psi_fit,psi_x*[1 1]);
    %}
end

if (Rx_pre-Rx)^2+(Zx_pre-Zx)^2>x_range^2/25
[Rx,Zx,psi_x,psi_fit] = fit_xpoint_initial(R,Z,psi,psi_lcs,psi_pri,Rx,Zx,x_range,icheck);
end

psi_fit=psi;%在初始化步骤中不准备改变初始的磁通分布
end

function [psi_x,Rx,Zx]=estimate_x(psi,psi_lcs,psi_pri,R,Z)
psi_pre=linspace(psi_lcs,psi_pri,25);
for irho=1:25%寻找x点的归一化磁通值
    h=contourc(R,Z,psi,psi_pre(irho)*[1,1]);
    if  h(2,1)==length(h)-1
        continue;
    elseif h(2,1)+h(2,h(2,1)+2)+2<length(h)%即等高线分为三段
        [m,I]=max(h(2,:));
        contour_core=h(:,(I+1):(I+m));
        [Zx,index]=min(contour_core(2,:));
        Rx=contour_core(1,index);
        %{
        figure
        plot(contour_core(1,:),contour_core(2,:))
        hold on
        scatter(Rx,Zx)
        %}
        psi_x=psi_pre(irho);
        %(Rx,Zx,npsix)为x点的坐标
        break
    end
end
end
