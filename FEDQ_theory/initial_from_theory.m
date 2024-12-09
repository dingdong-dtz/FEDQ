function eq=initial_from_theory(nr,nt)
icheck=0;
R1=5;
R2=7;
Z2=4.5;
A=-0.0005;
B=0.01;
rg=linspace(3.8,8,50*25);
zg=linspace(-7,5,150*25);
[R,Z]=meshgrid(rg,zg);
psi=theory_sample(R,Z);
Pprime0=-2./Z2.^2-8./R2.^2;
FFprime0=2*R1^2/Z2^2;

Rlim=[4 10 10 4 4];
Zlim=[-5 -5  2 2 -5];
R0=6.22;
Z0=-1.33;
psi_lcs=-0.8;
figure(1)
h=contour(rg,zg,psi,[0 0]+psi_lcs,'Visible','off');
rboundary=h(1,2:1+h(2,1));
zboundary=h(2,2:1+h(2,1));
psi_pri=-1.2;
h=contour(rg,zg,psi,[0 0]+psi_pri,'Visible','off');
rdown=h(1,2:1+h(2,1));
zdown=h(2,2:1+h(2,1));
close 1
[xup,yup,iiup] = polyxpoly(Rlim,Zlim,rboundary,zboundary);
[xdown,ydown,iidown] = polyxpoly(Rlim,Zlim,rdown,zdown);
divertor_left1=[xup(1) yup(1)];
divertor_left2=[xdown(1) ydown(1)];
divertor_right1=[xup(2) yup(2)];
divertor_right2=[xdown(2) ydown(2)];
divertor_left=divertor_left1+linspace(0,1,100)'*(divertor_left2-divertor_left1);
divertor_right=divertor_right1+linspace(0,1,100)'*(divertor_right2-divertor_right1);%记录顺序均是从上至下
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

rg=linspace(4,8,64*10);
zg=linspace(-7,5,128*10);

[R,Z]=meshgrid(rg,zg);
psi=theory_sample(R,Z)+0*rand(size(R));
%Initialize magnetic axis position
icheck=0;
axis_range=0.1;%Indicate the size of the area near the magnetic axis, regardless of accuracy during initialization
%[R0,Z0,dpsi_dr,psi_axis,psi_axis_fit]=findaxis(R,Z,rg,zg,psi,psi_lcs,icheck,R0,Z0,psi_axis);
[R0,Z0,dpsi_dr,psi_axis]=fit_axis_initial(R,Z,psi,R0,Z0,axis_range,icheck);
%Initialize x point position and magnetic flux value
icheck=0;
%[Rx,Zx,psi_x,psi] = fit_xpoint(R,Z,psi,psi_lcs,psi_pri,rg,zg,icheck);
Rx_pre=4.6;Zx_pre=-3.9;x_range=0.2;
[Rx,Zx,psi_x,psi] = fit_xpoint_initial(R,Z,psi,psi_lcs,psi_pri,Rx_pre,Zx_pre,x_range,icheck);

%Adjust the coordinate area range
Zmax=max(Zlcs);        [~,Zmax_index]=min(abs(zg-Zmax));
Zmin=min(Zpri);        [~,Zmin_index]=min(abs(zg-Zmin));
Rmax=max(Rlcs);        [~,Rmax_index]=min(abs(rg-Rmax));
Rmin=min([Rpri Rlcs]); [~,Rmin_index]=min(abs(rg-Rmin));
index1=max(1,Zmin_index-3):min(Zmax_index+3,length(zg));
index2=max(1,Rmin_index-3):min(Rmax_index+3,length(rg));
psi=psi(index1,index2);
R=R(index1,index2);
Z=Z(index1,index2);
%delta_psi=delta_psi(index1,index2);
rg=rg(index2);
zg=zg(index1);
%%
%Normalize the magnetic flux value
psib=psi_lcs-psi_axis;
npsi=(psi-psi_axis)/psib;
ndpsi_dr0=dpsi_dr/psib;
nrho=real(sqrt(npsi));%Prevent complex problems caused by fitting magnetic axis flux greater than the minimum flux value
npsi_x=(psi_x-psi_axis)/psib;               nrho_x=sqrt(npsi_x);
npsi_pri=(psi_pri-psi_axis)/psib;           nrho_pri=sqrt(npsi_pri);
npsi_left=(psi_left-psi_axis)/psib;         nrho_left=sqrt(npsi_left);  nrho_left(1) =1;    nrho_left(end) =nrho_pri;
npsi_right=(psi_right-psi_axis)/psib;       nrho_right=sqrt(npsi_right);nrho_right(1)=1;    nrho_right(end)=nrho_pri;
%%
%Rotational coordinate system
[alpha,Length]=cart2pol(R0-Rx,Z0-Zx);
[x,y] = rotate(R,Z,R0,Z0,alpha,Length);
[xdiv_left,ydiv_left]=rotate(Rdiv_left,Zdiv_left,R0,Z0,alpha,Length);
[xdiv_right,ydiv_right]=rotate(Rdiv_right,Zdiv_right,R0,Z0,alpha,Length);

%%
%Construct polar and radial coordinates
if ~exist("nr","var")
nr=64*4;
nt=64*4*4;
end
rho=x_nonuniform_new([0 1],nr,[0 nrho_x],[0.1 0.025],[3 20],0);
rho = x_adjustment(rho,[0 1],[nrho_pri,nrho_x]);
nr_down =sum(rho<nrho_pri)+1;
nr_inner=sum(rho<nrho_x);
nr_outer=nr-nr_inner;

theta=x_nonuniform_new([-pi/2 pi/2*5],nt,[0 2*pi 0 2*pi 4.9/4*pi 0 2*pi],[pi/24 pi/24 pi/48/4 pi/48/4 pi/6 pi/9 pi/9],[5 5 160*2 160*2 2 2 2]/1.5,0);
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
ng_r=3*max((t_br2-t_br1),1);
ng_l=3*max((t_bl1-t_bl2),1);
ng=5;% Representing the complete number of θ direction points participating in the reassignment of coordinates
%ng=min([10,nt1+1-t_br1-5-3,t_bl1-(nt1+nt_inner++1+5)-3]);% Representing the complete number of θ direction points participating in the reassignment of coordinates
nelse=0;% Indicate the number of additional nodes added
tr_min=t_br1-nelse;  tr_max=t_br2+ng_r;   nf_r=tr_max-tr_min; tf_r=(1:nf_r)/nf_r;
tl_max=t_bl1+nelse;  tl_min=t_bl2-ng_l;   nf_l=tl_max-tl_min; tf_l=(1:nf_l)/nf_l;
t_min=tr_min;   t_max=tl_max;%New polar coordinate boundary situation
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
    't_boundary',t_boundary,'tr_max',tr_max,'tl_min',tl_min);%,'r_boundary',r_boundary);

%%
paras.rho_inner=0.02;
paras.n_inner=sum(rho<paras.rho_inner);
%Initialize the magnetic surface coordinate system grid
x_new=zeros(nr,nt);
y_new=zeros(nr,nt);
R_new=zeros(nr,nt);
Z_new=zeros(nr,nt);
r_new=zeros(nr,nt);
%R_new(1,nt1+1:nt1+nt_inner)=R0;
%Z_new(1,nt1+1:nt1+nt_inner)=Z0;
R_new(1,:)=R0;
Z_new(1,:)=Z0;
figure(1)%The contour function requires a graph window
for j=2:nr
    h=contour(x,y,npsi,rho(j)^2*[1 1],'visible','off');%find surface positions
    %if j<paras.n_inner
    if j<=6
    r_new(j,nt1+1:nt1+nt_inner)=sqrt(rho(j)^2./(ndpsi_dr0(1)+(ndpsi_dr0(2)-ndpsi_dr0(1))*sin(theta_core+alpha-pi).^2));
    R_new(j,nt1+1:nt1+nt_inner)=r_new(j,nt1+1:nt1+nt_inner).*cos(theta_core+alpha-pi)+R0;
    Z_new(j,nt1+1:nt1+nt_inner)=r_new(j,nt1+1:nt1+nt_inner).*sin(theta_core+alpha-pi)+Z0;
    elseif j<=nr_down-1
        if h(2,1)<length(h)-1%There is not only one contour line obtained
            if h(2,2)<-1%The starting position is recorded below
                h(:,1:h(2,1)+1)=[];%Clear the useless section below
            end
        end
        jx1=h(1,2:h(2,1));%Due to being closed, the beginning and end of contour vectors are the same
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
        if max(h(2,2:h(2,1)+1))>0%(jx1, jy1) corresponds to the magnetic surface of the core, while (jx2, jy2) corresponds to the magnetic surface of the private area
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
        if j==nr_down%When j=nr_rown, the corresponding boundary line of the private area does not need to be solved again
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
            jr1=jr1(index_theta);
            r1=interp1([0 jtheta1 2*pi],[1 jr1 1],theta_core);
            [x1,y1] =  xy_theta_r_new(theta_core,r1);

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
            if max(h(2,2:h(2,1)+1))>0
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
    elseif j<=nr%For the SOL region, it is not possible to use θ, r coordinates for interpolation as it is too inaccurate. Instead, θ, x coordinates should be used for interpolation
        theta_f=theta;
        theta_f(tr_max-1:-1:tr_min)=tf_r*(theta_div_right_new(j)-theta(tr_min))+theta_f(tr_max-1:-1:tr_min);%The order of these points is from near the y-axis to both sides
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
close 1
[Theta,Rho]=meshgrid(theta,rho);
npsi=Rho.^2;
psi=npsi*psib+psi_axis;

%%
Pprime=zeros(size(rho))+Pprime0;
FFprime=zeros(size(rho))+FFprime0;
Pprime_pri=zeros(size(rho))+Pprime0;
FFprime_pri=zeros(size(rho))+FFprime0;
gpsi=rho.^2*psib+psi_axis;
gnpsi=rho.^2;
gpsi_pri=gpsi;
%%
%Set other parameters
paras.imax=40;
paras.err=1e-15;
%All parameters controlling the range of point X must be linked to the grid density, otherwise it is impossible to improve the calculation accuracy at point X during grid refinement
paras.X_y_range=(y_new(nr_inner+1-3,nt1+1)+1)*Length;
paras.X_x_range=(x_new(nr_inner+1+3,nt1+1))*Length;%Used to determine the range that needs to be updated with fitting when updating the grid
paras.xn_range=2;%Modify the xn_range surfaces on both sides of the interface using fitting methods
paras.theta_range=(theta(nt1+1+3)-theta(nt1+1));%Used to determine the weight of the fitted magnetic surface when updating the grid, making the magnetic surface at point X more accurate
paras.rho_range=rho(nr_inner+1)-rho(nr_inner+1-3);
eq = struct('R',R_new,'Z',Z_new,'theta',theta,'rho',rho,'geometry',geometry,'boundary',boundary, ...
    'R0',R0,'Z0',Z0,'Rx',Rx,'Zx',Zx,'paras',paras, ...
    'psi',psi,'psi_lcs',psi_lcs,'psib',psib,'psi_axis',psi_axis,'psi_pri',psi_pri,'psi_x',psi_x,...
    'ndpsi_dr',ndpsi_dr0,'npsi',npsi,'npsi_x',npsi_x,'npsi_left',npsi_left,'npsi_right',npsi_right,'npsi_pri',npsi_pri, ...
    'gpsi',gpsi,'gnpsi',gnpsi,'Pprime',Pprime,'FFprime',FFprime,'gpsi_pri',gpsi_pri,'Pprime_pri',Pprime_pri,'FFprime_pri',FFprime_pri);


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


