function r = r_xy(x,y,geometry)
%r = r_xy(x,y,geometry)
%   通过x和y来计算r，但是不需要计算θ，因为θ网格是选定的
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
r=zeros(nr,nt);
r(1:nr,nt1+1:nt1+nt_inner)=sqrt(x(1:nr,nt1+1:nt1+nt_inner).^2+(y(1:nr,nt1+1:nt1+nt_inner)).^2);
r(nr_down:nr,[t_min:nt1 nt1+nt_inner+1:t_max])=...
    sqrt(x(nr_down:nr,[t_min:nt1 nt1+nt_inner+1:t_max]).^2+(y(nr_down:nr,[t_min:nt1 nt1+nt_inner+1:t_max])+2).^2);

end