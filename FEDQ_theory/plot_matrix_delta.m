function plot_matrix_delta(R,Z,M,M_name,geometry,boundary,levels,with_color)
%M_name = inputname(3)%该矩阵的变量名称
nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nt=geometry.nt;
nr_inner=geometry.nr_inner;

nr_down=geometry.nr_down;
nr=geometry.nr;
t_min=geometry.t_min;
t_max=geometry.t_max;
if exist("levels","var")
    M_pre=levels;
else
    index=false(nr,nt);
    index(1:nr_down-1,nt1+1:nt1+nt_inner)=true;
    index(nr_down:nr,t_min:t_max)=true;
    Mmin=min(M(index));
    Mmax=max(M(index));
    if ~isfinite(Mmax)||~isfinite(Mmin)
        M_pre=25;
    else
        M_pre=linspace(Mmin,Mmax,25);
    end
end
%[R,Z,M] = change_for_contour(R,Z,M,geometry)
if exist('M_name','var')&&~isempty(M_name)
figure('Name',M_name,'NumberTitle','off');
end
hold on
if exist('with_color','var')&&with_color==0
    contour(R(1:nr_inner+1,[nt1+1:nt1+nt_inner,nt1+1]),Z(1:nr_inner+1,[nt1+1:nt1+nt_inner,nt1+1]),M(1:nr_inner+1,[nt1+1:nt1+nt_inner,nt1+1]),M_pre,"r-");
    contour(R(nr_down:nr_inner+1,[t_min:nt1,nt1+nt_inner+1:t_max]),Z(nr_down:nr_inner+1,[t_min:nt1,nt1+nt_inner+1:t_max]),M(nr_down:nr_inner+1,[t_min:nt1,nt1+nt_inner+1:t_max]),M_pre,"r-");
    contour(R(nr_inner+1:nr,t_min:t_max),Z(nr_inner+1:nr,t_min:t_max),M(nr_inner+1:nr,t_min:t_max),M_pre,"r-");
else
    contourf(R(1:nr_inner+1,[nt1+1:nt1+nt_inner,nt1+1]),Z(1:nr_inner+1,[nt1+1:nt1+nt_inner,nt1+1]),M(1:nr_inner+1,[nt1+1:nt1+nt_inner,nt1+1]),M_pre,'EdgeAlpha',0);
    contourf(R(nr_down:nr_inner+1,[t_min:nt1,nt1+nt_inner+1:t_max]),Z(nr_down:nr_inner+1,[t_min:nt1,nt1+nt_inner+1:t_max]),M(nr_down:nr_inner+1,[t_min:nt1,nt1+nt_inner+1:t_max]),M_pre,'EdgeAlpha',0);
    contourf(R(nr_inner+1:nr,t_min:t_max),Z(nr_inner+1:nr,t_min:t_max),M(nr_inner+1:nr,t_min:t_max),M_pre,'EdgeAlpha',0);
    colorbar
end
%%
%显示x点和过x点的几条网格线，包扩十字网格和等高线
scatter(R(1,nt1+1),Z(1,nt1+1),12,'red','filled');
scatter(R(nr_inner+1,nt1+1),Z(nr_inner+1,nt1+1),12,'red','filled');
%plot(R(:,nt1+1),Z(:,nt1+1),'r-')
%plot(R(nr_down:end,nt1+nt_inner+1),Z(nr_down:end,nt1+nt_inner+1),'r-')
plot(R(nr_inner+1,t_min:t_max),Z(nr_inner+1,t_min:t_max),'r-')

%%
if exist('boundary','var')&&isfield(boundary,'Rdiv_left')
plot(boundary.Rdiv_left,boundary.Zdiv_left,'k-');
plot(boundary.Rdiv_right,boundary.Zdiv_right,'k-');
plot(boundary.Rlcs,boundary.Zlcs,'k-');
plot(boundary.Rpri,boundary.Zpri,'k-');
end
axis equal
xlabel('R');
ylabel('Z');
box on
end