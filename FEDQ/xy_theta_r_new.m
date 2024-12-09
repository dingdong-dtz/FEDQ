function [x,y] = xy_theta_r_new(theta,radia,up_down)
%[x,y] = xy_theta_r(theta,radia)
%求出[theta,r]对应的坐标值（x，y）
%本程序的输出和输入是相同大小的矩阵，值一一对应
%up_down=1时表示核心区域的值，up_down=-1时表示私有区域的值，对于x点和SOL区域不需要该参数
if ~exist("up_down","var")
up_down=1;
end
x=zeros(size(radia));y=zeros(size(radia));
eqn=@(chi,r,tc) tan(pi/2*(-r*cos(chi)))+pi/2*(r*sin(chi))*cot(tc);
%x=r*cos(chi-pi/2)
%y=r*sin(chi-pi/2)
for n=1:length(theta(:))
    t=theta(n);
    r=radia(n);
    if isnan(r)%如果r为非数值量，即实际上的超出边界的量
        x_now=NaN;
        y_now=NaN;
    elseif t<0
        tc=-t;
        eq=@(chi) eqn(chi,r,tc);
        if r<1
            chi_s=fzero(eq,[tc,pi/2]);
        else
            chi_s=fzero(eq,[max(tc,acos(1/r)),pi/2]);
        end
        x_now=r*sin(chi_s);
        y_now=-r*cos(chi_s);
        y_now=-2-y_now;
    elseif t>0&&t<pi/2
        tc=t;
        eq=@(chi) eqn(chi,r,tc);
        if r<1
            chi_s=fzero(eq,[tc,pi/2]);
        else
            chi_s=fzero(eq,[max(tc,acos(1/r)),pi/2]);
        end
        x_now=r*sin(chi_s);
        y_now=-r*cos(chi_s);
    elseif t>pi/2*3&&t<pi*2
        tc=2*pi-t;
        eq=@(chi) eqn(chi,r,tc);
        if r<1
            chi_s=fzero(eq,[tc,pi/2]);
        else
            chi_s=fzero(eq,[max(tc,acos(1/r)),pi/2]);
        end
        chi_s=2*pi-chi_s;
        x_now=r*sin(chi_s);
        y_now=-r*cos(chi_s);
    elseif t>2*pi
        tc=t-2*pi;
        eq=@(chi) eqn(chi,r,tc);
        if r<1
            chi_s=fzero(eq,[tc,pi/2]);
        else
            chi_s=fzero(eq,[max(tc,acos(1/r)),pi/2]);
        end
        x_now=r*sin(chi_s);
        y_now=-r*cos(chi_s);
        x_now=-x_now;
        y_now=-2-y_now;
    elseif t>=pi/2&&t<=pi/2*3
        x_now=r*cos(t-pi/2);
        y_now=r*sin(t-pi/2);
    elseif t==0
        if r<1
            y_now=-1+up_down*(1-r);x_now=0;
        elseif r>1
            y_now=-1;x_now=sqrt(r^2-1);
        else
            y_now=-1;x_now=0;
        end
    elseif t==2*pi
        if r<1
            y_now=-1+up_down*(1-r);x_now=0;
        elseif r>1
            y_now=-1;x_now=-sqrt(r^2-1);
        else
            y_now=-1;x_now=0;
        end
    end
    x(n)=x_now;
    y(n)=y_now;
end
end