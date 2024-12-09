function Mmax = max_geometry(M,geometry)
%UNTITLED11 此处显示有关此函数的摘要
%   此处显示详细说明
nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nt=geometry.nt;
nr_down=geometry.nr_down;
nr=geometry.nr;
t_min=geometry.t_min;
t_max=geometry.t_max;
index=false(nr,nt);
index(1:nr_down-1,nt1+1:nt1+nt_inner)=true;
index(nr_down:nr,t_min:t_max)=true;
%Mmin=min(M(index));
Mmax=max(M(index));
end