function plot_contour_comparing(R,Z,Rx,Zx,geometry,boundary,psi_pre,psi_X_theory)
%plot_contour_comparing(R,Z,Rx,Zx,geometry,boundary,psi_pre)
%撰写文章需要的后处理程序
%本程序用于绘制划分的网格，并描绘边界
%本程序用于比较求解出的等磁通面和真实的等磁通面之间的差别
nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nt2=geometry.nt2;
nt=geometry.nt;
nr_inner=geometry.nr_inner;
nr_down=geometry.nr_down;
nr=geometry.nr;
t_min=geometry.t_min;
t_max=geometry.t_max;
Rmax=max_geometry(R,geometry);
Rmin=-max_geometry(-R,geometry);
Zmax=max_geometry(Z,geometry);
Zmin=-max_geometry(-Z,geometry);
lll=floor((nr-nr_inner)/3);
sss=mod(nr_inner+1,lll);
choose=lll:lll:nr;
choose=choose-sss;
index=find(choose==nr_inner+1);
choosenew=choose;
choosenew(index)=[];
fff=figure('Name','grid and boundary','NumberTitle','off');
hold on
f = @(x,y) theory_sample(x,y);

Linwide=0.5;
fcontour(f,[Rmin,Rmax,Zmin,Zmax],'LevelList',psi_pre(choosenew),'MeshDensity',500,'LineColor',"r","LineStyle","-","LineWidth",Linwide);
fcontour(f,[Rmin,Rmax,Zmin,Zmax],'LevelList',psi_X_theory,'MeshDensity',500,'LineColor',"r","LineStyle","-","LineWidth",Linwide*1.5);


%scatter(Rx,Zx,'blue','x','LineWidth',1);
scatter(R(1,nt1+1),Z(1,nt1+1),'blue','+','LineWidth',1);
j=3;
f1=plot(R(j,[nt1+1:nt1+nt_inner,nt1+1]),Z(j,[nt1+1:nt1+nt_inner,nt1+1]),'r-',"LineWidth",1);
f2=plot(R(j,[nt1+1:nt1+nt_inner,nt1+1]),Z(j,[nt1+1:nt1+nt_inner,nt1+1]),'b--',"LineWidth",1);
%这两个语句是为了产生图例而添加的
for j=choosenew%1:nr
    if j<=nr_inner&&j>=5
        plot(R(j,[nt1+1:nt1+nt_inner,nt1+1]),Z(j,[nt1+1:nt1+nt_inner,nt1+1]),'b--',"LineWidth",Linwide);
    end
    if j>=nr_down&&j<=nr_inner
        plot(R(j,[t_min:nt1 nt1+nt_inner+1:t_max]),Z(j,[t_min:nt1 nt1+nt_inner+1:t_max]),'b--',"LineWidth",Linwide);
    end
    if j>=nr_inner+1&&j<=nr
        plot(R(j,t_min:t_max),Z(j,t_min:t_max),'b--',"LineWidth",Linwide);
    end
end

plot(R(nr_inner+1,t_min:t_max),Z(nr_inner+1,t_min:t_max),'b--',"LineWidth",Linwide*1.5);
if exist('boundary','var')
    plot(boundary.Rdiv_left,boundary.Zdiv_left,'k-');
    plot(boundary.Rdiv_right,boundary.Zdiv_right,'k-');
    plot(boundary.Rlcs,boundary.Zlcs,'k-');
    plot(boundary.Rpri,boundary.Zpri,'k-');
end
fill([boundary.Rpri flip(boundary.Rpri)],[boundary.Zpri ones(size(boundary.Zpri))*(-5)],'w','EdgeAlpha',0);
axis equal
kuan=max(max(boundary.Rlcs),max(boundary.Rpri))-min(min(boundary.Rlcs),min(boundary.Rpri));
gao=-min(min(boundary.Zlcs),min(boundary.Zpri))+max(max(boundary.Zlcs),max(boundary.Zpri));
xlim([min(min(boundary.Rlcs),min(boundary.Rpri))-kuan*0.05 max(max(boundary.Rlcs),max(boundary.Rpri))+kuan*0.05])
ylim([min(min(boundary.Zlcs),min(boundary.Zpri))-gao*0.05 max(max(boundary.Zlcs),max(boundary.Zpri))+gao*0.05])
xlabel('R');
ylabel('Z');
box on
legend([f1,f2],'Analytical','Simulated','Location','best','FontSize',8)
%set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 12 9])
fff.Color=[1,1,1];
fff.Units="centimeters";
%fff.Position=[12,2,12,18];
%exportgraphics(gca,'compare300.png','Resolution',300)
kuanfig=4;
gaofig=kuanfig/kuan*gao;
set(gcf,'unit','centimeters','Position',[12,2,kuanfig+1+0.5,gaofig+1.2+0.5]);  
set(gca,'unit','centimeters','Position',[1 1.2 kuanfig gaofig]); 
exportgraphics(gcf,'compare.png','Resolution',1000)
%如上语句做出的图像时宽度为8cm的图像
end
