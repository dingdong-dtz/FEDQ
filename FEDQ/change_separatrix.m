%This program is used to draw images that indicate the application prospects of the program in a paper
%This program can input a series of eq and draw their boundaries on a graph to compare differences
%Using eq0 as a reference, draw the boundary conditions based on it
eq=init_from_grid;
eq1 = change_boundary(eq,1,+0.15);
eq2 = change_boundary(eq,1,-0.12);
eq0=FDEQ(eq);
eq1=FDEQ(eq1);
eq2=FDEQ(eq2);
R=eq0.R;
Z=eq0.Z;
geometry=eq0.geometry;
boundary=eq0.boundary;
nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nt=geometry.nt;
nr_inner=geometry.nr_inner;

nr_down=geometry.nr_down;
nr=geometry.nr;
t_min=geometry.t_min;
t_max=geometry.t_max;
%Display the x-point and several grid lines passing through the x-point, including the expanded cross grid and contour lines
fff=figure
hold on
scatter(R(1,nt1+1),Z(1,nt1+1));
f1=plot(R(nr_inner+1,t_min:t_max),Z(nr_inner+1,t_min:t_max),'r-',LineWidth=1.5);
%%
if exist('boundary','var')&&isfield(boundary,'Rdiv_left')
plot(boundary.Rdiv_left,boundary.Zdiv_left,'k-',LineWidth=1.5);
plot(boundary.Rdiv_right,boundary.Zdiv_right,'k-',LineWidth=1.5);
plot(boundary.Rlcs,boundary.Zlcs,'k-',LineWidth=1.5);
plot(boundary.Rpri,boundary.Zpri,'k-',LineWidth=1.5);
end
%%
Rmin=min(min(boundary.Rlcs),min(boundary.Rpri));
Rmax=max(max(boundary.Rlcs),max(boundary.Rpri));
Zmin=min(boundary.Zpri);
Zmax=max(boundary.Zlcs);
kuan=Rmax-Rmin;
gao=Zmax-Zmin;
%%
R=eq1.R;
Z=eq1.Z;
geometry=eq1.geometry;
boundary=eq1.boundary;
nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nt=geometry.nt;
nr_inner=geometry.nr_inner;

nr_down=geometry.nr_down;
nr=geometry.nr;
t_min=geometry.t_min;
t_max=geometry.t_max;
%Draw the changed boundary line
f2=plot(R(nr_inner+1,t_min:t_max),Z(nr_inner+1,t_min:t_max),'b-',LineWidth=1.5);
%%
R=eq2.R;
Z=eq2.Z;
geometry=eq2.geometry;
boundary=eq2.boundary;
nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nt=geometry.nt;
nr_inner=geometry.nr_inner;

nr_down=geometry.nr_down;
nr=geometry.nr;
t_min=geometry.t_min;
t_max=geometry.t_max;
f3=plot(R(nr_inner+1,t_min:t_max),Z(nr_inner+1,t_min:t_max),'g-',LineWidth=1.5);
%%
axis equal
xlabel('R(m)');
ylabel('Z(m)');
box on
xlim([Rmin-kuan*0.05,Rmax+kuan*0.05]);
ylim([Zmin-gao*0.05,Zmax+gao*0.05]);
legend([f1 f2 f3],{'separatrix-original ','separatrix-down','separatrix-up'},"FontSize",9)
fff.Color=[1,1,1];

kuanfig=8;
gaofig=kuanfig/kuan*gao;
set(gcf,'unit','centimeters','Position',[12,2,kuanfig+1.5+0.5,gaofig+1.2+0.5]);  
set(gca,'unit','centimeters','Position',[1.5 1.2 kuanfig gaofig]); 
exportgraphics(gca,'compare_with_EFIT.png','Resolution',1000)