function [x,y,npsi] = change_for_contour(x,y,npsi,geometry)
%[x,y,npsi] = change_for_contour(x,y,npsi,geometry)
%   Add appropriate partitions to the matrix so that contour lines can be represented as normal
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

index=false(nr,nt);
index(1:nr_down-1,nt1+1:nt1+nt_inner)=true;
index(nr_down:nr,t_min:t_max)=true;
npsi(~index)=NaN;

x_supply1   =[x(1:nr_inner,nt1+nt_inner+1);x(nr_inner+1:end,nt1)];
y_supply1   =[y(1:nr_inner,nt1+nt_inner+1);y(nr_inner+1:end,nt1)];
npsi_supply1=[npsi(1:nr_inner,nt1+nt_inner+1);npsi(nr_inner+1:end,nt1)];
x_supply2=[x(1:nr_inner,nt1+1);x(nr_inner+1:end,nt1+nt_inner)];
y_supply2=[y(1:nr_inner,nt1+1);y(nr_inner+1:end,nt1+nt_inner)];
npsi_supply2=[npsi(1:nr_inner,nt1+1);npsi(nr_inner+1:end,nt1+nt_inner)];
x_NaN1=x(:,nt1+1);
y_NaN1=y(:,nt1+1);
npsi_NaN1=npsi(:,nt1+1);
npsi_NaN1(1:nr_inner,1)=NaN;
x_NaN2=x(:,nt1+nt_inner+1);
y_NaN2=y(:,nt1+nt_inner+1);
npsi_NaN2=npsi(:,nt1+nt_inner+1);
npsi_NaN2(1:nr_inner,1)=NaN;
x=[x(:,t_min:nt1) x_supply1 x_NaN1  x(:,nt1+1:nt1+nt_inner) x_supply2 x_NaN2 x(:,nt1+nt_inner+1:t_max)];
y=[y(:,t_min:nt1) y_supply1 y_NaN1  y(:,nt1+1:nt1+nt_inner) y_supply2 y_NaN2 y(:,nt1+nt_inner+1:t_max)];
npsi=[npsi(:,t_min:nt1) npsi_supply1 npsi_NaN1  npsi(:,nt1+1:nt1+nt_inner) npsi_supply2 npsi_NaN2 npsi(:,nt1+nt_inner+1:t_max)];

end