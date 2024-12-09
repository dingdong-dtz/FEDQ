function plot_contour(R,Z,Rx,Zx,geometry,boundary)
%plot_contour(R,Z,Rx,Zx,geometry,boundary)
%本程序用于绘制划分的网格，并描绘边界
nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nt2=geometry.nt2;
nt=geometry.nt;
nr_inner=geometry.nr_inner;
nr_down=geometry.nr_down;
nr=geometry.nr;
t_min=geometry.t_min;
t_max=geometry.t_max;
figure('Name','grid and boundary','NumberTitle','off');
hold on
scatter(Rx,Zx);
scatter(R(1,nt1+1),Z(1,nt1+1));

for j=1:nr_inner
    plot(R(j,[nt1+1:nt1+nt_inner,nt1+1]),Z(j,[nt1+1:nt1+nt_inner,nt1+1]),'b.-')
end
for j=nr_down:nr_inner
    plot(R(j,[t_min:nt1 nt1+nt_inner+1:t_max]),Z(j,[t_min:nt1 nt1+nt_inner+1:t_max]),'r.-')
end
for j=nr_inner+1:nr
    plot(R(j,t_min:t_max),Z(j,t_min:t_max),'g.-')
end

if exist('boundary','var')
plot(boundary.Rdiv_left,boundary.Zdiv_left,'k-');
plot(boundary.Rdiv_right,boundary.Zdiv_right,'k-');
plot(boundary.Rlcs,boundary.Zlcs,'k-');
plot(boundary.Rpri,boundary.Zpri,'k-');
end
axis equal
end