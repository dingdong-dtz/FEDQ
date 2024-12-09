function theta = change_theta(nt,x_interval,f_interval)
%   theta = change_theta(nt,x_interval,f_interval)
%   该函数用于依据度规修改极向网格密度
L=length(x_interval);
Max_f_interval=max(f_interval);
y_edge_left=((0:floor(L/10)-1)/floor(L/10)).^4*(-Max_f_interval+f_interval(1))+Max_f_interval;
y_edge_right=((1-(1:floor(L/10))/floor(L/10))).^4*(-Max_f_interval+f_interval(end))+Max_f_interval;
x_edge_left=(0:floor(L/10)-1)/floor(L/10)*(x_interval(1)--pi/2)+-pi/2;
x_edge_right=(1:floor(L/10))/floor(L/10)*(pi*5/2-x_interval(end))+x_interval(end);
f_interval_new=[y_edge_left f_interval y_edge_right];
x_interval_new=[x_edge_left x_interval x_edge_right];
x_interval_new(1)=-pi/2;x_interval_new(end)=pi*5/2;

theta=x_nonuniform_real([-pi/2,pi*5/2],nt,x_interval_new,f_interval_new,1);
theta=x_adjustment(theta,[-pi/2 pi/2*5],[0 2*pi]);

end