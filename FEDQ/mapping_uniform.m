function ss = mapping_uniform(ipre,jpre,geometry)
%s = mapping_uniform(i,j,geometry)
%   将irho，jtheta网格以二元函数mapping映射到一个单调增加的数
%   要求以相对齐整的边界存储
%   本函数要求输入的i和j具有相同的元素数,和mapping_uniform的区别就在于此时ij一一对应构成坐标
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
t_all=t_max-t_min+1;
[j,i]=meshgrid(jpre,ipre);
sn = numel(i);
ss=zeros(size(i));
for kn=1:sn
    irho=i(kn);
    jtheta=j(kn);
    if irho==1
        if jtheta>nt1&&jtheta<=nt1+nt_inner
            s=1;
        else
            s=0;
            %warning('exceeding the grid range when mapping')
        end
    elseif irho<nr_down
        if jtheta>nt1&&jtheta<=nt1+nt_inner
            s=jtheta-nt1+(irho-2)*nt_inner+1;
        else
            s=0;
            %warning('exceeding the grid range when mapping')
        end
    elseif irho<=nr_inner
        if jtheta>nt1&&jtheta<=nt1+nt_inner
            s=jtheta-nt1+(irho-2)*nt_inner+1;
        elseif jtheta<=nt1&&jtheta>=t_min
            s=nt_inner*(nr_inner-1)+1+nr_outer*t_all...
                +jtheta-t_min+1+(irho-nr_down)*(nt1-t_min+1);
        elseif jtheta>nt1+nt_inner&&jtheta<=t_max
            s=nt_inner*(nr_inner-1)+1+nr_outer*t_all...
                +(nr_inner-nr_down+1)*(nt1-t_min+1)...
                +jtheta-nt1-nt_inner+(irho-nr_down)*(t_max-nt_inner-nt1);
        else
            s=0;
            %warning('exceeding the grid range when mapping')
        end
    elseif irho<=nr
        if jtheta>=t_min&&jtheta<=t_max
            s=nt_inner*(nr_inner-1)+1+jtheta-t_min+1+(irho-1-nr_inner)*t_all;
        else
            s=0;
            %warning('exceeding the grid range when mapping');
        end
    else
        s=0;
        %warning('exceeding the grid range when mapping');
    end
    ss(kn)=s;
end
end