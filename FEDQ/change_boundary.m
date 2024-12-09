function eq = change_boundary(eq,out_in,pp)
%UNTITLED2 此处显示有关此函数的摘要
%   out_in=1时表示修改内偏滤器靶板上的分布
%   pp表示靶板处分布偏移程度
geometry=eq.geometry;
nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nt=geometry.nt;
nr_inner=geometry.nr_inner;

nr_down=geometry.nr_down;
nr=geometry.nr;
t_min=geometry.t_min;
t_max=geometry.t_max;
if out_in==1
    R2=eq.R(nr_down,t_max);
    R1=eq.R(nr,t_max);

    Z2=eq.Z(nr_down,t_max);%偏滤器靶板左侧上方的端点
    Z1=eq.Z(nr,t_max);
    len=norm([R1-R2,Z1-Z2]);
    v_left=[R2-R1,Z2-Z1]/len;%单位化的偏滤器方向向量,由芯部指向下方
    u_left=[-v_left(2),v_left(1)];
    lwl=abs((eq.R(nr,t_min)-R1)*u_left(1)+(eq.Z(nr,t_min)-Z1)*u_left(2));%表示从外偏滤器下端点到左偏滤器所在直线的距离
else
    R2=eq.R(nr_down,t_min);
    R1=eq.R(nr,t_min);

    Z2=eq.Z(nr_down,t_min);%偏滤器靶板左侧上方的端点
    Z1=eq.Z(nr,t_min);
    len=norm([R1-R2,Z1-Z2]);
    v_left=[R2-R1,Z2-Z1]/len;%单位化的偏滤器方向向量,由芯部指向下方
    u_left=[-v_left(2),v_left(1)];
    lwl=abs((eq.R(nr,t_max)-R1)*u_left(1)+(eq.Z(nr,t_max)-Z1)*u_left(2));
    %{
l1=(eq.R-R1)*v_left(1)+(eq.Z-Z1)*v_left(2);%格点到与偏滤器靶板垂直的直线的距离
l2=abs((eq.R-R1)*u_left(1)+(eq.Z-Z1)*u_left(2));%格点到偏滤器靶板的距离
lwl=l2(nr,t_min);%表示从外偏滤器下端点到左偏滤器所在直线的距离
degree=exp(-l2/lwl*4).*(l1/len).*(l1>0)*pp;
eq.R=eq.R+degree*v_left(1);
eq.Z=eq.Z+degree*v_left(2);
    %}
end
%{
[eq.R,eq.Z] = change(eq.R,eq.Z,R1,Z1,v_left,u_left,len,lwl,pp);
[eq.boundary.Rpri,eq.boundary.Zpri] = change(eq.boundary.Rpri,eq.boundary.Zpri,R1,Z1,v_left,u_left,len,lwl,pp);
[eq.boundary.Rdiv_left,eq.boundary.Zdiv_left] = change(eq.boundary.Rdiv_left,eq.boundary.Zdiv_left,R1,Z1,v_left,u_left,len,lwl,pp);
[eq.boundary.Rdiv_right,eq.boundary.Zdiv_right] = change(eq.boundary.Rdiv_right,eq.boundary.Zdiv_right,R1,Z1,v_left,u_left,len,lwl,pp);
[eq.boundary.Rlcs,eq.boundary.Zlcs] = change(eq.boundary.Rlcs,eq.boundary.Zlcs,R1,Z1,v_left,u_left,len,lwl,pp);
%}
arrow=[eq.Rx-eq.R0,eq.Zx-eq.Z0];
al=norm(arrow);
arrow=arrow/al;
X_point=[eq.Rx,eq.Zx];
axis=[eq.R0,eq.Z0];
[eq.R,eq.Z] = change2(eq.R,eq.Z,X_point,[R1,Z1],arrow,u_left,al,lwl,pp);
[eq.boundary.Rpri,eq.boundary.Zpri] = change2(eq.boundary.Rpri,eq.boundary.Zpri,X_point,[R1,Z1],arrow,u_left,al,lwl,pp);
[eq.boundary.Rdiv_left,eq.boundary.Zdiv_left] = change2(eq.boundary.Rdiv_left,eq.boundary.Zdiv_left,X_point,[R1,Z1],arrow,u_left,al,lwl,pp);
[eq.boundary.Rdiv_right,eq.boundary.Zdiv_right] = change2(eq.boundary.Rdiv_right,eq.boundary.Zdiv_right,X_point,[R1,Z1],arrow,u_left,al,lwl,pp);
[eq.boundary.Rlcs,eq.boundary.Zlcs] = change2(eq.boundary.Rlcs,eq.boundary.Zlcs,X_point,[R1,Z1],arrow,u_left,al,lwl,pp);

end

function [Rxin,Zxin] = change(R,Z,R1,Z1,v,u,len,lwl,pp)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
l1=(R-R1)*v(1)+(Z-Z1)*v(2);%格点到与偏滤器靶板垂直的直线的距离
l2=abs((R-R1)*u(1)+(Z-Z1)*u(2));%格点到偏滤器靶板的距离
degree=exp(-l2/lwl*3).*(l1/len).*(l1>0)*pp;
Rxin=R+degree*v(1);
Zxin=Z+degree*v(2);
end
function [Rxin,Zxin] = change2(R,Z,v0,u0,v,u,vl,ul,pp)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
l1=(R-v0(1))*v(1)+(Z-v0(2))*v(2);%格点到与偏滤器靶板垂直的直线的距离
l2=abs((R-u0(1))*u(1)+(Z-u0(2))*u(2));%格点到偏滤器靶板的距离
degree=exp(-l2/ul*2).*(l1/vl*2).*(l1>0)*pp;
Rxin=R+degree*v(1);
Zxin=Z+degree*v(2);
end