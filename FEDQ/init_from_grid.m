function eq=init_from_grid(nr,nt,rg,zg,psi,Rlim,Zlim)
%This function initializes from the magnetic flux values on a rectangular grid
clearvars
if ~exist('psi','var')
    [rg,zg,psi,Rlim,Zlim,R0,Z0,psi_axis,eq_right,Rx,Zx]=data;
end
[R,Z]=meshgrid(rg,zg);
[dpsidR,dpsidZ]=cald_2d(psi,rg,zg);
[dpsidR2,~]=cald_2d(dpsidR./R,rg,zg);
[~,dpsidZ2]=cald_2d(dpsidZ,rg,zg);
delta_psi=dpsidR2.*R+dpsidZ2;
%%
figure(1)
hold on
%Determine a suitable outermost magnetic surface
psi_lcs=-0.318+0.001;
h=contour(rg,zg,psi,[psi_lcs,psi_lcs]);
rboundary=h(1,2:1+h(2,1));
zboundary=h(2,2:1+h(2,1));

%Determine a suitable private area magnetic surface
psi_pri=-0.331;
h=contour(rg,zg,psi,[psi_pri,psi_pri]);
rdown=h(1,2:1+h(2,1));
zdown=h(2,2:1+h(2,1));

p = 0.9;
value = csaps((1:length(rboundary)),[rboundary;zboundary],p,1:0.1:length(rboundary));
rboundary=value(1,:);
zboundary=value(2,:);
value = csaps((1:length(rdown)),[rdown;zdown],p,1:0.1:length(rdown));
rdown=value(1,:);
zdown=value(2,:);

psi_pri=-0.331;
close 1
[xup,yup,iiup] = polyxpoly(Rlim,Zlim,rboundary,zboundary);
[xdown,ydown,iidown] = polyxpoly(Rlim,Zlim,rdown,zdown);
divertor_left1=[xup(1) yup(1)];
divertor_left2=[xdown(1) ydown(1)];
divertor_right1=[xup(2) yup(2)];
divertor_right2=[xdown(2) ydown(2)];
divertor_left=divertor_left1+linspace(0,1,100)'*(divertor_left2-divertor_left1);
divertor_right=divertor_right1+linspace(0,1,100)'*(divertor_right2-divertor_right1);%The recording order is from top to bottom
psi_left =interp2(R,Z,psi,divertor_left(:,1),divertor_left(:,2));   psi_left(1)=psi_lcs;    psi_left(end)=psi_pri;
psi_right=interp2(R,Z,psi,divertor_right(:,1),divertor_right(:,2)); psi_right(1)=psi_lcs;   psi_right(end)=psi_pri;
Rdiv_left=divertor_left(:,1)';
Zdiv_left=divertor_left(:,2)';
Rdiv_right=divertor_right(:,1)';
Zdiv_right=divertor_right(:,2)';
Rlcs=[xup(1) rboundary(iiup(1,2)+1:iiup(2,2)) xup(2)];
Zlcs=[yup(1) zboundary(iiup(1,2)+1:iiup(2,2)) yup(2)];
Rpri=[xdown(1) rdown(iidown(1,2)+1:iidown(2,2)) xdown(2)];
Zpri=[ydown(1) zdown(iidown(1,2)+1:iidown(2,2)) ydown(2)];
[Rlcs,Zlcs] = interp2_curve(Rlcs,Zlcs,'linear',1000);
[Rpri,Zpri] = interp2_curve(Rpri,Zpri,'linear',250);


boundary=struct('Rlcs',Rlcs,'Zlcs',Zlcs, ...
    'Rpri',Rpri,'Zpri',Zpri, ...
    'Rdiv_left',Rdiv_left,'Rdiv_right',Rdiv_right,'Zdiv_left',Zdiv_left,'Zdiv_right',Zdiv_right, ...
    'psi_left',psi_left,'psi_right',psi_right,'psi_lcs',psi_lcs,'psi_pri',psi_pri);


%Initialize magnetic axis position
icheck=0;
[R0,Z0,dpsi_dr,psi_axis,~]=findaxis(R,Z,rg,zg,psi,psi_lcs,icheck,R0,Z0,psi_axis);
%Initialize x point position and magnetic flux value
icheck=0;
[Rx,Zx,psi_x,psi] = fit_xpoint(R,Z,psi,psi_lcs,psi_pri,rg,zg,icheck,Rx,Zx);

%%
%Adjust the coordinate area range
Zmax=max(Zlcs);        [~,Zmax_index]=min(abs(zg-Zmax));
Zmin=min(Zpri);        [~,Zmin_index]=min(abs(zg-Zmin));
Rmax=max(Rlcs);        [~,Rmax_index]=min(abs(rg-Rmax));
Rmin=min([Rpri Rlcs]); [~,Rmin_index]=min(abs(rg-Rmin));
psi=psi(Zmin_index-3:Zmax_index+3,Rmin_index-3:Rmax_index+3);
R=R(Zmin_index-3:Zmax_index+3,Rmin_index-3:Rmax_index+3);
Z=Z(Zmin_index-3:Zmax_index+3,Rmin_index-3:Rmax_index+3);
delta_psi=delta_psi(Zmin_index-3:Zmax_index+3,Rmin_index-3:Rmax_index+3);
rg=rg(Rmin_index-3:Rmax_index+3);
zg=zg(Zmin_index-3:Zmax_index+3);
%%
%Normalize the magnetic flux value
psib=psi_lcs-psi_axis;
npsi=(psi-psi_axis)/psib;
ndpsi_dr0=(dpsi_dr-psi_axis)/psib;
nrho=real(sqrt(npsi));%Prevent complex problems caused by fitting magnetic axis flux greater than the minimum flux value
npsi_x=(psi_x-psi_axis)/psib;               nrho_x=sqrt(npsi_x);
npsi_pri=(psi_pri-psi_axis)/psib;           nrho_pri=sqrt(npsi_pri);
npsi_left=(psi_left-psi_axis)/psib;         nrho_left=sqrt(npsi_left);  nrho_left(1) =1;    nrho_left(end) =nrho_pri;
npsi_right=(psi_right-psi_axis)/psib;       nrho_right=sqrt(npsi_right);nrho_right(1)=1;    nrho_right(end)=nrho_pri;%用于修正数值计算误差
%%
%Rotating coordinate system
[alpha,Length]=cart2pol(R0-Rx,Z0-Zx);
[x,y] = rotate(R,Z,R0,Z0,alpha,Length);
[xdiv_left,ydiv_left]=rotate(Rdiv_left,Zdiv_left,R0,Z0,alpha,Length);
[xdiv_right,ydiv_right]=rotate(Rdiv_right,Zdiv_right,R0,Z0,alpha,Length);

%%
%Construct polar and radial coordinates
if ~exist("nr","var")
nr=64*2;
nt=128*4;
end

rho=x_nonuniform_new([0 1],nr,[0 nrho_x],[0.1 0.025],[3 20],0);
rho = x_adjustment(rho,[0 1],[nrho_pri,nrho_x]);
nr_down =sum(rho<nrho_pri)+1;
nr_inner=sum(rho<nrho_x);
nr_outer=nr-nr_inner;

theta=x_nonuniform_new([-pi/2 pi/2*5],nt,[0 2*pi 0 2*pi 4.9/4*pi 0 2*pi],[pi/24 pi/24 pi/48/4 pi/48/4 pi/9 pi/9 pi/9],[5 5 40 40 2 2 2],0);
theta=x_adjustment(theta,[-pi/2 pi/2*5],[0 2*pi]);
nt1=sum(theta<0);
nt2=sum(theta>2*pi)+1;
nt_inner=nt-nt1-nt2;
theta_core=theta(nt1+1:nt-nt2);
theta_pri_right=theta(1:nt1);
theta_pri_left=theta(nt-nt2+1:end);
theta_pri=[theta_pri_right theta_pri_left];

theta_div_left  = coordinate_construct(xdiv_left, ydiv_left);
theta_div_right = coordinate_construct(xdiv_right,ydiv_right);% The θ value corresponding to the point on the polarizer corresponds one-to-one with the point recorded at the position of the polarizer
theta_div_left_new =interp1(nrho_left, theta_div_left, rho(nr_down:end));
theta_div_right_new=interp1(nrho_right,theta_div_right, rho(nr_down:end));% Indicate the boundary value of θ on each magnetic surface limited by the polarizer
theta_div_left_new =[ones(1,nr_down-1)*pi/2*5  theta_div_left_new];
theta_div_right_new=[ones(1,nr_down-1)*(-pi/2) theta_div_right_new];% Convert it into a form with corresponding values for all magnetic surfaces to facilitate the following calculations

t_boundary=zeros(3,nr);
for j=1:nr_down-1
    t_boundary(1,j)=nt1+1;
    t_boundary(2,j)=nt1+nt_inner;
    t_boundary(3,j)=nt_inner;
end
for j=nr_down:nr
    t_boundary(1,j)=sum(theta<theta_div_right_new(j))+1;
    t_boundary(2,j)=sum(theta<theta_div_left_new(j));
    t_boundary(3,j)=t_boundary(2,j)-t_boundary(1,j)+1;
end
t_br1=min(t_boundary(1,nr_down:end));   t_br2=max(t_boundary(1,nr_down:end));
t_bl1=max(t_boundary(2,nr_down:end));   t_bl2=min(t_boundary(2,nr_down:end));
ng_r=min(3*max((t_br2-t_br1),1),10);
ng_l=min(3*max((t_bl1-t_bl2),1),10);
nelse=0;
tr_min=t_br1-nelse;  tr_max=t_br2+ng_r;   nf_r=tr_max-tr_min; tf_r=(1:nf_r)/nf_r;
tl_max=t_bl1+nelse;  tl_min=t_bl2-ng_l;   nf_l=tl_max-tl_min; tf_l=(1:nf_l)/nf_l;
t_min=tr_min;   t_max=tl_max;
t_boundary(1,1:nr_down-1)=nt1+1;        t_boundary(1,nr_down:end)=t_min;
t_boundary(2,1:nr_down-1)=nt1+nt_inner; t_boundary(2,nr_down:end)=t_max;
t_boundary(3,1:nr_down-1)=nt_inner;     t_boundary(3,nr_down:end)=t_max-t_min+1;
theta_pri_f_all=zeros(nr_inner+1-nr_down+1,t_max-nt_inner-t_min+1);

tf_r=tf_r.^2;
tf_l=tf_l.^2;
for j=nr_down:nr_inner+1
    theta_pri_f=theta_pri;
    theta_pri_f(tr_max-1:-1:tr_min)=tf_r*(theta_div_right_new(j)-theta_pri(tr_min))+theta_pri_f(tr_max-1:-1:tr_min);%这一段点的顺序由靠近y轴向两边
    theta_pri_f(tl_min-nt_inner+1:tl_max-nt_inner)=tf_l*(theta_div_left_new(j)-theta_pri(tl_max-nt_inner))+theta_pri_f(tl_min-nt_inner+1:tl_max-nt_inner);
    theta_pri_f_all(j-nr_down+1,:)=theta_pri_f(t_min:t_max-nt_inner);
end
geometry=struct('nr',nr,'nr_inner',nr_inner,'nr_outer',nr_outer,'nr_down',nr_down, ...
    'nt',nt,'nt_inner',nt_inner,'nt1',nt1,'nt2',nt2,'t_min',t_min,'t_max',t_max,'alpha',alpha,'Length',Length,...
    't_boundary',t_boundary,'tr_max',tr_max,'tl_min',tl_min);

%%
paras.rho_inner=0.02;
paras.n_inner=sum(rho<paras.rho_inner);
%close all
%Initialize the magnetic surface coordinate system grid
x_new=zeros(nr,nt);
y_new=zeros(nr,nt);
R_new=zeros(nr,nt);
Z_new=zeros(nr,nt);
r_new=zeros(nr,nt);
R_new(1,:)=R0;
Z_new(1,:)=Z0;
for j=2:nr
    h=contour(x,y,npsi,rho(j)^2*[1 1],'visible','off');%find surface positions
    if j<=9
    r_new(j,nt1+1:nt1+nt_inner)=sqrt(rho(j)^2./(ndpsi_dr0(1)+(ndpsi_dr0(2)-ndpsi_dr0(1))*sin(theta_core+alpha-pi).^2));
    R_new(j,nt1+1:nt1+nt_inner)=r_new(j,nt1+1:nt1+nt_inner).*cos(theta_core+alpha-pi)+R0;
    Z_new(j,nt1+1:nt1+nt_inner)=r_new(j,nt1+1:nt1+nt_inner).*sin(theta_core+alpha-pi)+Z0;
    elseif j<=nr_down-1
        if h(2,1)<length(h)-1
            if h(2,2)<-1
                h(:,1:h(2,1)+1)=[];
            end
        end
        jx1=h(1,2:h(2,1));
        jy1=h(2,2:h(2,1));
        jtheta1=coordinate_construct(jx1,jy1);
        
        [jtheta1,index]=sort(jtheta1);jx1=jx1(index);jy1=jy1(index);
        curve = csape([ jtheta1 jtheta1(1)+2*pi],[ jx1 jx1(1); jy1 jy1(1)],'period');
        value = fnval(curve,mod(theta_core-jtheta1(1),2*pi)+jtheta1(1));x1=value(1,:);y1=value(2,:);
        
        x_new(j,nt1+1:nt1+nt_inner)=x1;
        y_new(j,nt1+1:nt1+nt_inner)=y1;
        [R_new(j,nt1+1:nt1+nt_inner),Z_new(j,nt1+1:nt1+nt_inner)] = rotate_inverse(x1,y1,R0,Z0,alpha,Length);
    elseif j<=nr_inner
        theta_pri_f=theta_pri_f_all(j-nr_down+1,:);
        if max(h(2,2:h(2,1)+1))>0
            jx1=h(1,2:h(2,1)+1);
            jy1=h(2,2:h(2,1)+1);
            jx2=h(1,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
            jy2=h(2,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
        else
            jx2=h(1,2:h(2,1)+1);
            jy2=h(2,2:h(2,1)+1);
            jx1=h(1,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
            jy1=h(2,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
        end
        if j==nr_down
            [jx2,jy2]=rotate(Rpri,Zpri,R0,Z0,alpha,Length);
        end
        jtheta1=coordinate_construct(jx1,jy1);
        jtheta2=coordinate_construct(jx2,jy2);

        [jtheta1,index]=sort(jtheta1);jx1=jx1(index);jy1=jy1(index);
        curve = csape([ jtheta1 jtheta1(1)+2*pi],[ jx1 jx1(1); jy1 jy1(1)],'period');
        value = fnval(curve,mod(theta_core-jtheta1(1),2*pi)+jtheta1(1));x1=value(1,:);y1=value(2,:);
        curve = csape(mod(jtheta2+pi/2,2*pi)-pi/2,[ jx2;jy2]);
        value = fnval(curve,mod(theta_pri_f+pi,2*pi)-pi);x2=value(1,:);y2=value(2,:);

        x_new(j,nt1+1:nt1+nt_inner)=x1;
        y_new(j,nt1+1:nt1+nt_inner)=y1;
        x_new(j,[t_min:nt1,nt1+nt_inner+1:t_max])=x2;
        y_new(j,[t_min:nt1,nt1+nt_inner+1:t_max])=y2;
        [R_new(j,t_min:t_max),Z_new(j,t_min:t_max)] =...
            rotate_inverse(x_new(j,t_min:t_max),y_new(j,t_min:t_max),R0,Z0,alpha,Length);

    elseif j==nr_inner+1
        theta_pri_f=theta_pri_f_all(j-nr_down+1,:);
        if h(2,1)==length(h)-1
            jx=h(1,2:h(2,1)+1);
            jy=h(2,2:h(2,1)+1);
            para_pre=0.1:0.1:0.9;
            jx1=jx(jy>-1&abs(jx)>0.025);jx1=[para_pre*jx1(1) jx1 para_pre*jx1(end)];
            jy1=jy(jy>-1&abs(jx)>0.025);jy1=[para_pre*(jy1(1)+1)-1 jy1 para_pre*(jy1(end)+1)-1];
            jx2=jx(jy<-1&abs(jx)>0.05);[jx2,index_x]=sort(jx2);  number_left=sum(jx2<0);
            jy2=jy(jy<-1&abs(jx)>0.05);jy2=jy2(index_x);
            jx2=[jx2(1:number_left) para_pre*jx2(number_left)       0   para_pre*jx2(number_left+1)         jx2(number_left+1:end)];
            jy2=[jy2(1:number_left) para_pre*(jy2(number_left)+1)-1 -1  para_pre*(jy2(number_left+1)+1)-1   jy2(number_left+1:end)];

            jtheta1=coordinate_construct(jx1,jy1);
            [jtheta1,index_theta]=sort(jtheta1);jx1=jx1(index_theta);jy1=jy1(index_theta);
            jr1=sqrt(jx1.^2+jy1.^2);
            r1=interp1([0 jtheta1 2*pi],[1 jr1 1],theta_core);
            [x1,y1] =  xy_theta_r_new(theta_core,r1,1);

            jtheta2=coordinate_construct(jx2,jy2);jtheta2(number_left+length(para_pre)+1)=0;
            jr2=sqrt(jx2.^2+(jy2+2).^2);
            r2=interp1(mod(jtheta2+pi/2,2*pi)-pi/2,jr2,mod(theta_pri_f+pi/2,2*pi)-pi/2);
            [x2,y2] =  xy_theta_r_new(theta_pri_f,r2);

            x_new(j,nt1+1:nt1+nt_inner)=x1;
            y_new(j,nt1+1:nt1+nt_inner)=y1;
            x_new(j,[t_min:nt1,nt1+nt_inner+1:t_max])=x2;
            y_new(j,[t_min:nt1,nt1+nt_inner+1:t_max])=y2;
            [R_new(j,t_min:t_max),Z_new(j,t_min:t_max)] =...
                rotate_inverse(x_new(j,t_min:t_max),y_new(j,t_min:t_max),R0,Z0,alpha,Length);
        else
            if max(h(2,2:h(2,1)+1))>0%Number 1 corresponds to the magnetic surface of the core, and number 2 corresponds to the magnetic surface of the private area
                jx_up=h(1,2:h(2,1)+1);
                jy_up=h(2,2:h(2,1)+1);
                jx_down=h(1,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
                jy_down=h(2,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
            else
                jx_down=h(1,2:h(2,1)+1);
                jy_down=h(2,2:h(2,1)+1);
                jx_up=h(1,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
                jy_up=h(2,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
            end
            para_pre=0.1:0.1:0.9;
            jx1=jx_up(jy_up>-1&abs(jx_up)>0.025);jx1=[ para_pre*jx1(1) jx1 para_pre*jx1(end) ];
            jy1=jy_up(jy_up>-1&abs(jx_up)>0.025);jy1=[ para_pre*(jy1(1)+1)-1 jy1 para_pre*(jy1(end)+1)-1 ];
            jx2=jx_down(jy_down<-1&abs(jx_down)>0.05);[jx2,index_x]=sort(jx2);  number_left=sum(jx2<0);
            jy2=jy_down(jy_down<-1&abs(jx_down)>0.05);jy2=jy2(index_x);
            jx2=[jx2(1:number_left) para_pre*jx2(number_left)       0   para_pre*jx2(number_left+1)         jx2(number_left+1:end)];
            jy2=[jy2(1:number_left) para_pre*(jy2(number_left)+1)-1 -1  para_pre*(jy2(number_left+1)+1)-1   jy2(number_left+1:end)];

            jtheta1=coordinate_construct(jx1,jy1);[jtheta1,index_theta]=sort(jtheta1);jx1=jx1(index_theta);jy1=jy1(index_theta);
            jtheta2=coordinate_construct(jx2,jy2);jtheta2(number_left+length(para_pre)+1)=0;
            index_right=number_left+length(para_pre)+1:length(jtheta2);index_left=1:number_left+length(para_pre)+1;
            
            curve1 = csape([0 jtheta1 2*pi],[0 jx1 0;-1 jy1 -1]);
            value1 = fnval(curve1,mod(theta_core-jtheta1(1),2*pi)+jtheta1(1));x1=value1(1,:);y1=value1(2,:);
            index_right_new=1:nt1-t_min+1;index_left_new=nt1-t_min+2:t_max-t_min+1-nt_inner;
            curve2_right = csape(mod(jtheta2(index_right)+pi,2*pi)-pi,[jx2(index_right);jy2(index_right)]);
            curve2_left =  csape(mod(jtheta2(index_left) +pi,2*pi)-pi,[jx2(index_left );jy2(index_left )]);
            value2_right= fnval(curve2_right,mod(theta_pri_f(index_right_new)+pi,2*pi)-pi);x2_right=value2_right(1,:);y2_right=value2_right(2,:);
            value2_left = fnval(curve2_left ,mod(theta_pri_f(index_left_new )+pi,2*pi)-pi);x2_left= value2_left(1,:); y2_left =value2_left (2,:);
            x2=[x2_right x2_left];y2=[y2_right y2_left];

            x_new(j,nt1+1:nt1+nt_inner)=x1;                 x_new(j,nt1+1)=0;
            y_new(j,nt1+1:nt1+nt_inner)=y1;                 y_new(j,nt1+1)=-1;
            x_new(j,[t_min:nt1,nt1+nt_inner+1:t_max])=x2;   x_new(j,nt1+nt_inner+1)=0;
            y_new(j,[t_min:nt1,nt1+nt_inner+1:t_max])=y2;   y_new(j,nt1+nt_inner+1)=-1;

            [R_new(j,t_min:t_max),Z_new(j,t_min:t_max)] =...
                rotate_inverse(x_new(j,t_min:t_max),y_new(j,t_min:t_max),R0,Z0,alpha,Length);
        end
    elseif j<=nr
        theta_f=theta;
        theta_f(tr_max-1:-1:tr_min)=tf_r*(theta_div_right_new(j)-theta(tr_min))+theta_f(tr_max-1:-1:tr_min);%这一段点的顺序由靠近y轴向两边
        theta_f(tl_min+1:tl_max)=tf_l*(theta_div_left_new(j)-theta(tl_max))+theta_f(tl_min+1:tl_max);
        theta_f=theta_f(t_min:t_max);
        if j==nr
            [jx1,jy1]=rotate(Rlcs,Zlcs,R0,Z0,alpha,Length);
        else
            jx1=h(1,2:h(2,1)+1);
            jy1=h(2,2:h(2,1)+1);
        end
        jtheta1=coordinate_construct(jx1,jy1);
        curve = csape(jtheta1,[ jx1;jy1]);

        value = fnval(curve,theta_f);x1=value(1,:);y1=value(2,:);
        x_new(j,t_min:t_max)=x1;y_new(j,t_min:t_max)=y1;
        [R_new(j,t_min:t_max),Z_new(j,t_min:t_max)] =...
            rotate_inverse(x_new(j,t_min:t_max),y_new(j,t_min:t_max),R0,Z0,alpha,Length);
    else
    end
    if j<nr_down
        R_new(j,1:nt1+1-1)=R_new(j,nt1+1);  Z_new(j,1:nt1+1-1)=Z_new(j,nt1+1);
        R_new(j,nt1+nt_inner+1:end)=R_new(j,nt1+nt_inner);Z_new(j,nt1+nt_inner+1:end)=Z_new(j,nt1+nt_inner);
    else
        R_new(j,1:t_min-1)=R_new(j,t_min);  Z_new(j,1:t_min-1)=Z_new(j,t_min);
        R_new(j,t_max+1:end)=R_new(j,t_max);Z_new(j,t_max+1:end)=Z_new(j,t_max);
    end
end
[R_new,Z_new]=smooth_1Dxin(R_new,Z_new,rho,theta,geometry,sqrt(npsi_x),0);
[~,r_rho]=meshgrid(theta,rho);
npsi=r_rho.^2;
psi=npsi*psib+psi_axis;
nd=4;% (nd+1) is the number of lattice points used in the difference scheme for constructing the system of equations, and an even number is generally used to ensure central symmetry
md=9;% MD is the number of points used to calculate the metric difference format, usually using odd numbers to ensure center symmetry
order=4;% MD is the accuracy order of the differential format for calculating metric differences, usually using odd numbers to ensure center symmetry
X_n_fit=5;
M=Metric_SOL_steady(R_new,Z_new,rho,theta,geometry,md,order,X_n_fit,0);%求解度规
%%
%Calculate voltage and current terms based on initial equilibrium
Pprime_new=zeros(size(rho));
FFprime_new=zeros(size(rho));
Pprime_pri_new=zeros(size(rho));
FFprime_pri_new=zeros(size(rho));
delta_new=interp2(R,Z,delta_psi,R_new,Z_new);
for j=2:nr
    if j<nr_inner
        [p,~]=polyfit(R_new(j,nt1+1:nt1+nt_inner).^2,delta_new(j,nt1+1:nt1+nt_inner),1);
        Pprime_new(j)=-p(1);
        FFprime_new(j)=-p(2);
    end
    if j>=nr_inner
        [p,~]=polyfit(R_new(j,t_min:t_max).^2,delta_new(j,t_min:t_max),1);
        Pprime_new(j)=-p(1);
        FFprime_new(j)=-p(2);
    end
    if j<=nr_inner&&j>=nr_down
        [p,~]=polyfit(R_new(j,[t_min:nt1,nt1+nt_inner+1:t_max]).^2,delta_new(j,[t_min:nt1,nt1+nt_inner+1:t_max]),1);
        Pprime_pri_new(j)=-p(1);
        FFprime_pri_new(j)=-p(2);
    end
end
index=7:25;
Pprime_new(1:15)=polyval(polyfit(rho(index).^2,Pprime_new(index),1),rho(1:15).^2);
FFprime_new(1:15)=polyval(polyfit(rho(index).^2,FFprime_new(index),1),rho(1:15).^2);

smoothpara=0.1;
Pprime_new(end)=0;
FFprime_new(end)=0;
Pprime_new=smoothdata(Pprime_new,2,'sgolay',smoothpara,'SamplePoints',rho,'degree',2);
FFprime_new=smoothdata(FFprime_new,2,'sgolay',smoothpara,'SamplePoints',rho,'degree',2);
Pprime_new(end)=0;
FFprime_new(end)=0;
Pprime_new=smoothdata(Pprime_new,2,'sgolay',smoothpara,'SamplePoints',rho,'degree',2);
FFprime_new=smoothdata(FFprime_new,2,'sgolay',smoothpara,'SamplePoints',rho,'degree',2);
Pprime_new(end)=0;
FFprime_new(end)=0;
Pprime_new=smoothdata(Pprime_new,2,'sgolay',smoothpara,'SamplePoints',rho,'degree',2);
FFprime_new=smoothdata(FFprime_new,2,'sgolay',smoothpara,'SamplePoints',rho,'degree',2);
Pprime_new(end)=0;
FFprime_new(end)=0;

Pprime=Pprime_new;
FFprime=FFprime_new;
Pprime_pri=Pprime_pri_new*0;
FFprime_pri=FFprime_pri_new*0;
gpsi=rho.^2*psib+psi_axis;
gnpsi=rho.^2;
gpsi_pri=gpsi;
Js=evalJs(gpsi,Pprime,FFprime,gpsi,Pprime_pri,FFprime_pri,R_new,psi,geometry);

%%
%Set other parameters
paras.imax=40;
paras.err=1e-15;
icheck=0;
paras.x_range=0.05*Length;% The size of the quadratic approximation region near point x in RZ space
paras.xn_range=2;% Modify the xn_range surfaces on both sides of the interface using fitting methods

paras.X_y_range=(y_new(nr_inner+1-3,nt1+1)+1)*Length;
paras.X_x_range=(x_new(nr_inner+1+3,nt1+1))*Length;
paras.theta_range=(theta(nt1+1+3)-theta(nt1+1));
paras.rho_range=rho(nr_inner+1)-rho(nr_inner+1-3);


eq = struct('R',R_new,'Z',Z_new,'theta',theta,'rho',rho,'geometry',geometry,'boundary',boundary, ...
    'R0',R0,'Z0',Z0,'Rx',Rx,'Zx',Zx,'paras',paras, ...
    'psi',psi,'psi_lcs',psi_lcs,'psib',psib,'psi_axis',psi_axis,'psi_pri',psi_pri,'psi_x',psi_x,...
    'ndpsi_dr',ndpsi_dr0,'npsi',npsi,'npsi_x',npsi_x,'npsi_left',npsi_left,'npsi_right',npsi_right,'npsi_pri',npsi_pri, ...
    'M',M,'Js',Js,'gpsi',gpsi,'gnpsi',gnpsi,'Pprime',Pprime,'FFprime',FFprime,'gpsi_pri',gpsi_pri,'Pprime_pri',Pprime_pri,'FFprime_pri',FFprime_pri);
end


function [R0,Z0,dpsi_dr,psi_axis,psi_axis_fit]=findaxis(R,Z,rg,zg,psi,psi_lcs,icheck,R0,Z0,psi_axis)
%**************************************
% function [R0,Z0,npsi,psibm,psi_axis,sIp,r0,t0]=findaxis(eq,icheck)
% find magnetic axis
% psi_axis is the smaller one within 'psi_axis_fit' and the smallest number of all psi
%The equipotential surface at the magnetic axis is close to an ellipse, and the position of the magnetic axis is determined by fitting the inner magnetic surface
%*************************************
if ~exist('psi_axis','var')
    psi_axis_pre=min(psi,[],'all');
    [index_R0,index_Z0] = find(psi==psi_axis_pre);
    R0_pre=R(index_R0,index_Z0);
    Z0_pre=Z(index_R0,index_Z0);
else
    psi_axis_pre=psi_axis;
    R0_pre=R0;
    Z0_pre=Z0;
    [~,index_R0]=min(abs(R0-rg));
    [~,index_Z0]=min(abs(Z0-zg));
end
near_R=R(index_Z0+(-5:5),index_R0+(-5:5));
near_Z=Z(index_Z0+(-5:5),index_R0+(-5:5));
near_psi=psi(index_Z0+(-5:5),index_R0+(-5:5));
%{
figure
contourf(near_R,near_Z,near_psi);
%}
roundEqn = '(a*(x-R0).^2+b*(y-Z0).^2)+psi_axis';
startPoints = [R0_pre,Z0_pre,1,1,0];
f1 = fit([near_R(:),near_Z(:)],near_psi(:),roundEqn,'Start', startPoints);
%coeffnames(f1)
R0=f1.R0;
Z0=f1.Z0;
dpsi_dr=[f1.a f1.b];
psi_axis_fit=f1.psi_axis;
if psi_axis_fit<psi_axis_pre
    psi_axis=psi_axis_fit;
else
    psi_axis=psi_axis_pre;
end
end

function [Rx,Zx,psi_x,psi_fit] = fit_xpoint(R,Z,psi,psi_lcs,psi_pri,rg,zg,icheck,Rx_pre,Zx_pre)
%[parameter,beta] = fit_xpoint(eq,icheck)
%Fit x point to find its position and magnetic flux value
if ~exist("Rx_pre","var")
[psi_x_pre,Rx_pre,Zx_pre]=estimate_x(psi,psi_lcs,psi_pri,R,Z,rg,zg);
end
[~,index_Rx]=min(abs(rg-Rx_pre));
[~,index_Zx]=min(abs(zg-Zx_pre));

near_R=R(index_Zx+(-5:5),index_Rx+(-5:5));
near_Z=Z(index_Zx+(-5:5),index_Rx+(-5:5));
near_psi=psi(index_Zx+(-5:5),index_Rx+(-5:5));
lambda=(rg(index_Rx+5)-rg(index_Rx-5))/2;%The quantity used to measure the size of the fitting area
sf = fit([near_R(:), near_Z(:)],near_psi(:),'poly22');
M=[sf.p20,sf.p11/2;
    sf.p11/2,sf.p02];
T=[sf.p10 sf.p01];
[V,D]=eig(M);%M*V = V*D
e1=V(:,1)/(V(:,1)'*V(:,1));
e2=V(:,2)/(V(:,2)'*V(:,2));
H=[e1,e2];
TH=T*H;
psi_x=sf.p00-TH(1)^2/4/D(1,1)-TH(2)^2/4/D(2,2);
X_Point=H*[-TH(1)/2/D(1,1);-TH(2)/2/D(2,2)];
Rx=X_Point(1);
Zx=X_Point(2);
psi_fit_diff=(sf(R,Z)-psi).*exp(-((R-Rx).^2+(Z-Zx).^2)/lambda^2*4);
psi_fit=psi+psi_fit_diff;
if icheck
    near_psi_fit=sf(near_R,near_Z);
    figure('Name','psi_fit','NumberTitle','off')
    hold on
    contour(near_R,near_Z,near_psi);
    contour(near_R,near_Z,near_psi_fit);
    scatter(Rx_pre,Zx_pre);
    scatter(Rx,Zx);
    %{
    figure
    hold on
    contour(R,Z,psi,linspace(psi_pri,psi_lcs,25));
    contour(R,Z,psi_fit,linspace(psi_pri,psi_lcs,25));
    contour(R,Z,psi,psi_x*[1 1]);
    contour(R,Z,psi_fit,psi_x*[1 1]);
    %}
end
end
function [psi_x,Rx,Zx]=estimate_x(psi,psi_lcs,psi_pri,R,Z,rg,zg)
psi_pre=linspace(psi_lcs,psi_pri,25);
for irho=1:25
    h=contourc(rg,zg,psi,psi_pre(irho)*[1,1]);
    if  h(2,1)==length(h)-1
        continue;
    elseif h(2,1)+h(2,h(2,1)+2)+2<length(h)
        [m,I]=max(h(2,:));
        contour_core=h(:,(I+1):(I+m));
        [Zx,index]=min(contour_core(2,:));
        Rx=contour_core(1,index);
        %{
        figure
        plot(contour_core(1,:),contour_core(2,:))
        hold on
        scatter(Rx,Zx)
        %}
        psi_x=psi_pre(irho);
        break
    end
end
end

function [rg,zg,psi,Rlim,Zlim,R0,Z0,psi_axis,eq_right,Rx,Zx]=data
load('data.mat');
psi=paras.eqdata.psirz;
rg=paras.eqdata.rg;
zg=paras.eqdata.zg;
R0=paras.R0;
Z0=paras.Z0;
Rx=paras.Rx;
Zx=paras.Zx;
psi_axis=paras.eqdata.psiaxis;
Rlim=paras.eqdata.Rlim;
Zlim=paras.eqdata.Zlim;
Pprime=paras.Pprime;
FFprime=paras.FFprime;
gpsi=paras.gpsi;
psiedge=paras.eqdata.psiedge;
gnpsi=gpsi;
gpsi=gpsi*(psiedge-psi_axis)+psi_axis;
eq_right=struct('Pprime',Pprime,'FFprime',FFprime,'gpsi',gpsi,'gnpsi',gnpsi);
end

function [x,y] = rotate(R,Z,R0,Z0,alpha,L)
%[x,y] = rotate(R,Z,R0,Z0,alpha,L)
%Rotate the coordinate system clockwise by pi/2- α, with (R0, Z0) as the center and reduce it by L times
x=((R-R0)*sin(alpha)-(Z-Z0)*cos(alpha))/L;
y=((R-R0)*cos(alpha)+(Z-Z0)*sin(alpha))/L;
end
function [R,Z] = rotate_inverse(x,y,R0,Z0,alpha,L)
%[R,Z] = rotate_inverse(x,y,R0,Z0,alpha,L)
%Enlarge by L times around (0,0) and then rotate the coordinate system counterclockwise by pi/2- α
R=(cos(alpha)*y+sin(alpha)*x)*L+R0;
Z=(sin(alpha)*y-cos(alpha)*x)*L+Z0;
end



