function eq=upRZ_SOL_old(eq,icheck)


x_n_fit=5;
[eq.R0,eq.Z0,eq.psi_axis,M0,eq.npsi,eq.psib]=fit_axis_new(eq.R,eq.Z,eq.psi,eq.geometry,eq.rho);
[eq.Rx,eq.Zx,eq.npsi_x,eq.npsi,SV,M,T,r_range] = fit_xpoint_new(eq.Rx,eq.Zx,eq.R,eq.Z,eq.npsi,eq.geometry,icheck,x_n_fit);
eq.r_range=r_range;
r_range
psib=eq.psib;
npsi=eq.npsi;
psi_axis=eq.psi_axis;
R0=eq.R0;
Z0=eq.Z0;
Rx=eq.Rx;
Zx=eq.Zx;
npsi_x=eq.npsi_x;nrho_x=sqrt(npsi_x);
eq.psi_x=npsi_x*psib+psi_axis;
%

npsi_pri=(eq.boundary.psi_pri-psi_axis)/psib;           nrho_pri=sqrt(npsi_pri);
npsi_left=(eq.boundary.psi_left-psi_axis)/psib;         nrho_left=sqrt(npsi_left);  nrho_left(1) =1;    nrho_left(end) =nrho_pri;
npsi_right=(eq.boundary.psi_right-psi_axis)/psib;       nrho_right=sqrt(npsi_right);nrho_right(1)=1;    nrho_right(end)=nrho_pri;%用于修正数值计算误差

nt=eq.geometry.nt;
nt1=eq.geometry.nt1;
nt2=eq.geometry.nt2;
nt_inner=eq.geometry.nt_inner;

nr=eq.geometry.nr;
rho=x_nonuniform_new([0 1],nr,[0 nrho_x nrho_x],[0.1 1-nrho_x (1-nrho_x)/10],[3 20 40],0);
rho = x_adjustment(rho,[0 1],[nrho_pri,nrho_x]);
nr_down =sum(rho<nrho_pri)+1;
nr_inner=sum(rho<nrho_x);
nr_outer=nr-nr_inner;


theta=eq.theta;
theta_core=theta(nt1+1:nt1+nt_inner);
theta_pri=[theta(1:nt1),theta(nt1+nt_inner+1:end)];
[alpha,Length]=cart2pol(R0-Rx,Z0-Zx);
[x,y] = rotate(eq.R,eq.Z,R0,Z0,alpha,Length);
[xdiv_left,ydiv_left]=rotate(eq.boundary.Rdiv_left,eq.boundary.Zdiv_left,R0,Z0,alpha,Length);
[xdiv_right,ydiv_right]=rotate(eq.boundary.Rdiv_right,eq.boundary.Zdiv_right,R0,Z0,alpha,Length);
theta_div_left  = coordinate_construct(xdiv_left, ydiv_left);
theta_div_right = coordinate_construct(xdiv_right,ydiv_right);% The θ value corresponding to the point on the polarizer corresponds one-to-one with the point recorded at the position of the polarizer
theta_div_left_new =interp1(nrho_left, theta_div_left, rho(nr_down:end),'pchip','extrap');
theta_div_right_new=interp1(nrho_right,theta_div_right, rho(nr_down:end),'pchip','extrap');% Indicate the boundary value of θ on each magnetic surface limited by the polarizer
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
nelse=0;% Indicate the number of additional nodes added
tr_min=t_br1-nelse;   tr_max=t_br2+ng_r;    nf_r=tr_max-tr_min;
tf_r=(1:nf_r)/nf_r;% Eliminate jump points by changing θ at the edge
tl_max=t_bl1+nelse;   tl_min=t_bl2-ng_l;    nf_l=tl_max-tl_min;
tf_l=(1:nf_l)/nf_l;% Eliminate jump points by changing θ at the edge
t_min=tr_min;    t_max=tl_max;% New polar coordinate boundary situation
t_boundary(1,1:nr_down-1)=nt1+1;        t_boundary(1,nr_down:end)=t_min;
t_boundary(2,1:nr_down-1)=nt1+nt_inner; t_boundary(2,nr_down:end)=t_max;
t_boundary(3,1:nr_down-1)=nt_inner;     t_boundary(3,nr_down:end)=t_max-t_min+1;
theta_pri_f_all=zeros(nr_inner+1-nr_down+1,t_max-nt_inner-t_min+1);
tf_r=tf_r.^2;
tf_l=tf_l.^2;
for j=nr_down:nr_inner+1
    theta_pri_f=theta_pri;
    theta_pri_f(tr_max-1:-1:tr_min)=tf_r*(theta_div_right_new(j)-theta_pri(tr_min))+theta_pri_f(tr_max-1:-1:tr_min);
    theta_pri_f(tl_min-nt_inner+1:tl_max-nt_inner)=tf_l*(theta_div_left_new(j)-theta_pri(tl_max-nt_inner))+theta_pri_f(tl_min-nt_inner+1:tl_max-nt_inner);
    theta_pri_f_all(j-nr_down+1,:)=theta_pri_f(t_min:t_max-nt_inner);
end

U=Length*[sin(alpha) cos(alpha)
    -cos(alpha) sin(alpha)];
Mxin=U'*M*U;
Txin=([2*Rx 2*Zx]*M+T)*U*0;
S=[0            0           Mxin(1,1)
    0            2*Mxin(1,2) Txin(1)
    Mxin(2,2)    Txin(2)     npsi_x];
M0xin=U'*M0*U;

theta_range=eq.paras.theta_range;
rho_range=eq.paras.rho_range;
rho_x=sqrt(npsi_x);
pp=1;
ww=1000;
[Theta,Rho]=meshgrid(theta,rho);
%%
nrho=real(sqrt(npsi));% Prevent complex problems caused by fitting magnetic axis flux greater than the minimum flux value
[x, y,nrho] = change_for_contour(x,y,nrho,eq.geometry);% The stored matrix segmentation is converted into a format that can calculate contour lines normally
x_new=zeros(nr,nt);
y_new=zeros(nr,nt);
R_new=zeros(nr,nt);
Z_new=zeros(nr,nt);
r_new=zeros(nr,nt);
X_x_range=r_range/Length;
X_y_range=r_range/Length;
x_n_range=3;


flux_rho=rho;%This quantity represents the selected value of the magnetic surface function used to search for the magnetic surface
for j=1:nr
    h=contour(x,y,nrho,rho(j)*[1 1],'visible','off');%find surface positions
    if j==1
        R_new(1,nt1+1:nt1+nt_inner)=R0;
        Z_new(1,nt1+1:nt1+nt_inner)=Z0;
    elseif j<5%The innermost contour line is fitted using elliptical contour lines
        r_new(j,nt1+1:nt1+nt_inner)=sqrt(rho(j)^2./(M0xin(1,1)*sin(theta_core).^2+M0xin(2,2)*cos(theta_core).^2-2*M0xin(1,2)*sin(theta_core).*cos(theta_core)));
        x_new(j,nt1+1:nt1+nt_inner)=r_new(j,nt1+1:nt1+nt_inner).*cos(theta_core-pi/2);
        y_new(j,nt1+1:nt1+nt_inner)=r_new(j,nt1+1:nt1+nt_inner).*sin(theta_core-pi/2);
        [R_new(j,nt1+1:nt1+nt_inner),Z_new(j,nt1+1:nt1+nt_inner)] = rotate_inverse(x_new(j,nt1+1:nt1+nt_inner),y_new(j,nt1+1:nt1+nt_inner),R0,Z0,alpha,Length);
    elseif j<=nr_down-1
        if sum(h(1,:)==flux_rho(j))>1
            h(:,h(1,:)==flux_rho(j))=[];
            jx1=h(1,h(2,:)>-1);
            jy1=h(2,h(2,:)>-1);%(jx1, jy1) corresponds to the magnetic surface of the core
        else
            jx1=h(1,2:h(2,1));%Due to being closed, the beginning and end of contour vectors are the same
            jy1=h(2,2:h(2,1));
        end
        jtheta1=coordinate_construct(jx1,jy1);
        [jtheta1,index_j]=sort(jtheta1);jx1=jx1(index_j);jy1=jy1(index_j);%%Reordering is for better interpolation
        jtheta1=[jtheta1(end-20:end)-2*pi jtheta1 jtheta1(1:20)+2*pi];

        weight=ones(size(jtheta1));
        weight=weight+ww*exp(-(rho(j)-rho_x).^2/rho_range^2).*exp(-(mod(jtheta1+pi,2*pi)-pi).^2/theta_range^2);
        p=ones(size(jtheta1))*pp;
        p=p+(1-pp)*exp(-(rho(j)-rho_x).^2/rho_range^2).*exp(-(mod(jtheta1+pi,2*pi)-pi).^2/theta_range^2);
        value=csaps(jtheta1,[jx1(end-20:end) jx1 jx1(1:20);jy1(end-20:end) jy1 jy1(1:20)],...
            p,theta_core,weight);
        x1=value(1,:);y1=value(2,:);
        x_new(j,nt1+1:nt1+nt_inner)=x1;y_new(j,nt1+1:nt1+nt_inner)=y1;
        [R_new(j,nt1+1:nt1+nt_inner),Z_new(j,nt1+1:nt1+nt_inner)] = rotate_inverse(x1,y1,R0,Z0,alpha,Length);
    elseif j<=nr_inner
        theta_pri_f=theta_pri_f_all(j-nr_down+1,:);
        if sum(h(1,:)==flux_rho(j))>2
            h(:,h(1,:)==flux_rho(j))=[];
            jx1=h(1,h(2,:)>-1);
            jy1=h(2,h(2,:)>-1);%(jx1, jy1) corresponds to the magnetic surface of the core, and (jx2, jy2) corresponds to the magnetic surface of the private area
            jx2=h(1,h(2,:)<-1);[jx2,index_j]=sort(jx2,"descend");
            jy2=h(2,h(2,:)<-1);jy2=jy2(index_j);
        elseif sum(h(1,:)==flux_rho(j))==2
            if max(h(2,2:h(2,1)+1))>0
                jx1=h(1,2:h(2,1)+1);
                jy1=h(2,2:h(2,1)+1);
                jx2=h(1,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));   [jx2,index_j]=sort(jx2,"descend");
                jy2=h(2,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));   jy2=jy2(index_j);
            else
                jx2=h(1,2:h(2,1)+1);    [jx2,index_j]=sort(jx2,"descend");
                jy2=h(2,2:h(2,1)+1);    jy2=jy2(index_j);
                jx1=h(1,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
                jy1=h(2,h(2,1)+3:h(2,1)+2+h(2,h(2,1)+2));
            end
        else
            jx1=h(1,2:h(2,1)+1);
            jy1=h(2,2:h(2,1)+1);
            [jx2,jy2]=rotate(eq.boundary.Rpri,eq.boundary.Zpri,R0,Z0,alpha,Length);
            [jx2,index2]=sort(jx2,'descend');jy2=jy2(index2);
            if j>nr_down  disp('error contour');  end
        end
        [theta1,~]=cart2pol(jx1,jy1);[~,index1]=sort(mod(theta1+pi/2,2*pi));
        jx1=jx1(index1);jy1=jy1(index1);
        [jx1,jy1]= remove_repeat(jx1,jy1);
        [jx2,jy2]= remove_repeat(jx2,jy2);

        if j==nr_down
            [jx2,jy2]=rotate(eq.boundary.Rpri,eq.boundary.Zpri,R0,Z0,alpha,Length);
            [jx2,index2]=sort(jx2,'descend');jy2=jy2(index2);
        elseif j>=nr_inner+1-x_n_range
            jx_else=linspace(-X_x_range,X_x_range,50);
            jy1_else=interp1(jx1(jy1<-0.5),jy1(jy1<-0.5),jx_else,'pchip','extrap');
            jy2_else=interp1(jx2,jy2,jx_else,'pchip','extrap');
            jx1=[jx1 jx_else];jy1=[jy1 jy1_else];
            jx2=[jx2 jx_else];jy2=[jy2 jy2_else];[jx2,index2]=sort(jx2,'descend');jy2=jy2(index2);
            index1=(abs(jx1)<X_x_range*2&jy1<0);
            index2=(abs(jx2)<X_x_range*2);
            P1=[jx1(index1)'.^2 jx1(index1)' ones(size(jx1(index1)'))]*S;
            %P(1)*x^2+P(2)*x+P(3)=psi
            delta1=sqrt(P1(:,2).^2-4*P1(:,1).*(P1(:,3)-rho(j)^2));
            jy1_else_fit=(-P1(:,2)+sign(Mxin(2,2))*delta1)./P1(:,1)/2-1;
            P2=[jx2(index2)'.^2 jx2(index2)' ones(size(jx2(index2)'))]*S;
            delta2=sqrt(P2(:,2).^2-4*P2(:,1).*(P2(:,3)-rho(j)^2));
            jy2_else_fit=(-P2(:,2)-sign(Mxin(2,2))*delta2)./P2(:,1)/2-1;

            jy1(index1)=jy1(index1)+(jy1_else_fit'-jy1(index1)).*exp(-2*(jx1(index1).^2/X_x_range^2+(jy1(index1)+1).^2/X_y_range^2).^2);
            jy2(index2)=jy2(index2)+(jy2_else_fit'-jy2(index2)).*exp(-2*(jx2(index2).^2/X_x_range^2+(jy2(index2)+1).^2/X_y_range^2).^2);
        end
        jtheta1=coordinate_construct(jx1,jy1);
        [jtheta1,index_j]=sort(jtheta1);jx1=jx1(index_j);jy1=jy1(index_j);
        [jtheta1,index_repeat]= remove_repeat1(jtheta1);jx1(index_repeat)=[];jy1(index_repeat)=[];
        jtheta2=coordinate_construct(jx2,jy2);
        jtheta1=[jtheta1(end-20:end)-2*pi jtheta1 jtheta1(1:20)+2*pi];

        weight=ones(size(jtheta1));
        weight=weight+ww*exp(-(rho(j)-rho_x).^2/rho_range^2).*exp(-(mod(jtheta1+pi,2*pi)-pi).^2/theta_range^2);
        p=ones(size(jtheta1))*pp;
        p=p+(1-pp)*exp(-(rho(j)-rho_x).^2/rho_range^2).*exp(-(mod(jtheta1+pi,2*pi)-pi).^2/theta_range^2);
        value=csaps(jtheta1,[jx1(end-20:end) jx1 jx1(1:20);jy1(end-20:end) jy1 jy1(1:20)],...
            p,theta_core,weight);
        x1=value(1,:);y1=value(2,:);
        jtheta2=mod(jtheta2+pi,2*pi)-pi;

        weight=ones(size(jtheta2));weight(1)=ww;weight(end)=ww;
        weight=weight+ww*exp(-(rho(j)-rho_x).^2/rho_range^2).*exp(-(mod(jtheta2+pi,2*pi)-pi).^2/theta_range^2);
        p=ones(size(jtheta2))*pp;
        p=p+(1-pp)*exp(-(rho(j)-rho_x).^2/rho_range^2).*exp(-(mod(jtheta2+pi,2*pi)-pi).^2/theta_range^2);
        value=csaps(jtheta2,[jx2;jy2],p,mod(theta_pri_f+pi/2,2*pi)-pi/2,weight);
        x2=value(1,:);y2=value(2,:);
        x_new(j,nt1+1:nt1+nt_inner)=x1;                 y_new(j,nt1+1:nt1+nt_inner)=y1;
        x_new(j,[t_min:nt1 nt1+nt_inner+1:t_max])=x2;   y_new(j,[t_min:nt1 nt1+nt_inner+1:t_max])=y2;

        x_new(j,t_min)=jx2(1,1);        y_new(j,t_min)=jy2(1,1);
        x_new(j,t_max)=jx2(1,end);      y_new(j,t_max)=jy2(1,end);
        [R_new(j,t_min:t_max),Z_new(j,t_min:t_max)] =...
            rotate_inverse(x_new(j,t_min:t_max),y_new(j,t_min:t_max),R0,Z0,alpha,Length);
    elseif j==nr_inner+1%Used to solve contour lines containing x points
        theta_pri_f=theta_pri_f_all(j-nr_down+1,:);
        SV=SV*[sin(alpha) cos(alpha);-cos(alpha) sin(alpha)];
        if SV(1,1)*SV(1,2)>0
            k1=-SV(2,1)/SV(2,2);%The lines representing the top right and bottom left
            k2=-SV(1,1)/SV(1,2);
        else
            k2=-SV(2,1)/SV(2,2);
            k1=-SV(1,1)/SV(1,2);
        end
        index_begin=h(1,:)==flux_rho(j);
        h(:,index_begin)=[];
        jx=h(1,:);jy=h(2,:);
        [jx,jy]= remove_repeat(jx,jy);
        para_pre=0.05:0.05:0.95;
        jx_up=jx(jy>-1&~(abs(jx)<X_x_range/sqrt(2)&jy<0));% (abs (jx)<X_x_range/2&jy<0) indicates the point near point x that needs to be deleted
        jy_up=jy(jy>-1&~(abs(jx)<X_x_range/sqrt(2)&jy<0));% This is to remove the incorrect points found on the climbing surface near point x, and to select this area to encrypt the coordinate points
        [theta_j,~]=cart2pol(jx_up,jy_up); theta_j=mod(theta_j+pi/2,2*pi);
        [~,index_j]=sort(theta_j); jx_up=jx_up(index_j); jy_up=jy_up(index_j);% Sort by θ coordinate again

        jx_up=[para_pre*jx_up(1) jx_up para_pre*jx_up(end)];
        jy_up=[para_pre*(jy_up(1)+1)-1 jy_up para_pre*(jy_up(end)+1)-1];
        jy_up_fit=k1*jx_up.*(jx_up>=0)+k2*jx_up.*(jx_up<0)-1;

        jy_up=jy_up+(jy_up_fit-jy_up).*exp(-2*(jx_up.^2/X_x_range^2+(jy_up+1).^2/X_y_range^2).^2);
        jtheta1=coordinate_construct(jx_up,jy_up);
        [jtheta1,index_theta]=sort(jtheta1);jx_up=jx_up(index_theta);jy_up=jy_up(index_theta);
        [jtheta1,index_repeat]= remove_repeat1(jtheta1);jx_up(index_repeat)=[];jy_up(index_repeat)=[];

        jx_down=jx(jy<-1&abs(jx)>X_x_range/sqrt(2));[jx_down,index_x]=sort(jx_down);  number_left=sum(jx_down<0);
        jy_down=jy(jy<-1&abs(jx)>X_x_range/sqrt(2));jy_down=jy_down(index_x);
        jx_down=[jx_down(1:number_left) para_pre*jx_down(number_left)       0   para_pre*jx_down(number_left+1)         jx_down(number_left+1:end)];
        jy_down=[jy_down(1:number_left) para_pre*(jy_down(number_left)+1)-1 -1  para_pre*(jy_down(number_left+1)+1)-1   jy_down(number_left+1:end)];
        jy_down_fit=k2*jx_down.*(jx_down>=0)+k1*jx_down.*(jx_down<0)-1;

        jy_down=jy_down+(jy_down_fit-jy_down).*exp(-2*(jx_down.^2/X_x_range^2+(jy_down+1).^2/X_y_range^2).^2);
        jtheta2=coordinate_construct(jx_down,jy_down);jtheta2(number_left+length(para_pre)+1)=0;

        jtheta1=[0 jtheta1 2*pi];
        weight=ones(size(jtheta1));
        weight=weight+ww*exp(-(rho(j)-rho_x).^2/rho_range^2).*exp(-(mod(jtheta1+pi,2*pi)-pi).^2/theta_range^2);
        p=ones(size(jtheta1))*pp;
        p=p+(1-pp)*exp(-(rho(j)-rho_x).^2/rho_range^2).*exp(-(mod(jtheta1+pi,2*pi)-pi).^2/theta_range^2);
        value=csaps(jtheta1,[0 jx_up 0;-1 jy_up -1],p,theta_core,weight);
        x1=value(1,:);y1=value(2,:);
        jtheta2=mod(jtheta2+pi,2*pi)-pi;

        weight=ones(size(jtheta2));
        weight=weight+ww*exp(-(rho(j)-rho_x).^2/rho_range^2).*exp(-(mod(jtheta2+pi,2*pi)-pi).^2/theta_range^2);
        weight=ones(size(jtheta2));weight(1)=ww;weight(end)=ww;
        p=ones(size(jtheta2))*pp;
        p=p+(1-pp)*exp(-(rho(j)-rho_x).^2/rho_range^2).*exp(-(mod(jtheta2+pi,2*pi)-pi).^2/theta_range^2);
        value=csaps(jtheta2,[jx_down;jy_down],1,mod(theta_pri_f+pi/2,2*pi)-pi/2,weight);
        x2=value(1,:);y2=value(2,:);
        x_new(j,nt1+1:nt1+nt_inner)=x1;                 y_new(j,nt1+1:nt1+nt_inner)=y1;
        x_new(j,[t_min:nt1,nt1+nt_inner+1:t_max])=x2;   y_new(j,[t_min:nt1,nt1+nt_inner+1:t_max])=y2;

        x_new(j,t_min)=jx_down(1,end);   y_new(j,t_min)=jy_down(1,end);
        x_new(j,t_max)=jx_down(1,1);   y_new(j,t_max)=jy_down(1,1);

        [R_new(j,t_min:t_max),Z_new(j,t_min:t_max)] =...
            rotate_inverse(x_new(j,t_min:t_max),y_new(j,t_min:t_max),R0,Z0,alpha,Length);
    elseif j<=nr
        if j==nr
            [jx,jy]=rotate(eq.boundary.Rlcs,eq.boundary.Zlcs,R0,Z0,alpha,Length);
            jtheta=coordinate_construct(jx,jy);
            [~,index]=sort(jtheta);
            jx=jx(index);jy=jy(index);
        else
            index_begin=h(1,:)==flux_rho(j);
            h(:,index_begin)=[];
            jx=h(1,:);jy=h(2,:);
            [jx,jy]= remove_repeat(jx,jy);
            if j<=nr_inner+1+x_n_range%Fit contour lines at point x
                jy_else=linspace(-X_y_range,X_y_range,50)-1;
                jx1_else=interp1(jy(jy<-0.5&jx<0),jx(jy<-0.5&jx<0),jy_else,'pchip','extrap');
                jx2_else=interp1(jy(jy<-0.5&jx>0),jx(jy<-0.5&jx>0),jy_else,'pchip','extrap');
                jx=[jx jx1_else jx2_else];jy=[jy jy_else jy_else];
                jtheta=coordinate_construct(jx,jy);
                [~,index]=sort(jtheta);
                jx=jx(index);jy=jy(index);
                index1=(jx<0&((jy+1).^2/X_y_range^2+jx.^2/X_x_range^2)<5^2);
                index2=(jx>0&((jy+1).^2/X_y_range^2+jx.^2/X_x_range^2)<5^2);
                P1=[(jy(index1)'+1).^2 jy(index1)'+1 ones(size(jy(index1)'))]*S';
                %P(1)*x^2+P(2)*x+P(3)=psi
                delta1=sqrt(P1(:,2).^2-4*P1(:,1).*(P1(:,3)-rho(j)^2));
                jx1_else_fit=(-P1(:,2)+sign(Mxin(2,2))*delta1)./P1(:,1)/2;
                P2=[(jy(index2)'+1).^2 jy(index2)'+1 ones(size(jy(index2)'))]*S';
                delta2=sqrt(P2(:,2).^2-4*P2(:,1).*(P2(:,3)-rho(j)^2));
                jx2_else_fit=(-P2(:,2)-sign(Mxin(2,2))*delta2)./P2(:,1)/2;

                jx(index1)=jx(index1)+(jx1_else_fit'-jx(index1)).*exp(-2*((jy(index1)+1).^2/X_y_range^2+jx(index1).^2/X_x_range^2).^2);
                jx(index2)=jx(index2)+(jx2_else_fit'-jx(index2)).*exp(-2*((jy(index2)+1).^2/X_y_range^2+jx(index2).^2/X_x_range^2).^2);

            end
        end
        jtheta1=coordinate_construct(jx,jy);
        [jtheta1,index]=sort(jtheta1);
        jx=jx(index);jy=jy(index);
        theta_f=theta;
        theta_f(tr_max-1:-1:tr_min)=tf_r*(theta_div_right_new(j)-theta(tr_min))+theta_f(tr_max-1:-1:tr_min);%这一段点的顺序由靠近y轴向两边
        theta_f(tl_min+1:tl_max)=tf_l*(theta_div_left_new(j)-theta(tl_max))+theta_f(tl_min+1:tl_max);

        weight=ones(size(jtheta1));
        weight=weight+ww*exp(-(rho(j)-rho_x).^2/rho_range^2).*exp(-(mod(jtheta1+pi,2*pi)-pi).^2/theta_range^2);
        p=ones(size(jtheta1))*pp;
        p=p+(1-pp)*exp(-(rho(j)-rho_x).^2/rho_range^2).*exp(-(mod(jtheta1+pi,2*pi)-pi).^2/theta_range^2);
        value=csaps(jtheta1,[jx;jy],p,theta_f(t_min:t_max),weight);

        jr=sqrt(jx.^2+min(jy.^2,(jy+2).^2));
        jrxin=csaps(jtheta1,jr,p,theta_f(t_min:t_max),weight);
        [x2,y2]=xy_theta_r_new(theta_f(t_min:t_max),jrxin);

        x2xin=csaps(jtheta1,jx,p,theta_f(t_min:t_max),weight);
        index=theta_f(t_min:t_max)<pi/2|theta_f(t_min:t_max)>pi/2*3;

        y2xin=y_theta_x(theta_f(t_min:t_max),x2xin);
        x2(index)=x2xin(index);
        y2(index)=y2xin(index);
        x_new(j,t_min:t_max)=x2;y_new(j,t_min:t_max)=y2;
        if jtheta1(end)>jtheta1(1)
            x_new(j,t_min)=jx(1,1);     y_new(j,t_min)=jy(1,1);
            x_new(j,t_max)=jx(1,end);   y_new(j,t_max)=jy(1,end);
        else
            x_new(j,t_min)=jx(1,end);   y_new(j,t_min)=jy(1,end);
            x_new(j,t_max)=jx(1,1);     y_new(j,t_max)=jy(1,1);
        end
        [R_new(j,t_min:t_max),Z_new(j,t_min:t_max)] =...
            rotate_inverse(x_new(j,t_min:t_max),y_new(j,t_min:t_max),R0,Z0,alpha,Length);
    end
    if j<nr_down
        R_new(j,1:nt1)=R_new(j,nt1+1);  Z_new(j,1:nt1)=Z_new(j,nt1+1);
        R_new(j,nt1+nt_inner+1:end)=R_new(j,nt1+nt_inner);Z_new(j,nt1+nt_inner+1:end)=Z_new(j,nt1+nt_inner);
    else
        R_new(j,1:t_min-1)=R_new(j,t_min);  Z_new(j,1:t_min-1)=Z_new(j,t_min);
        R_new(j,t_max+1:end)=R_new(j,t_max);Z_new(j,t_max+1:end)=Z_new(j,t_max);
    end
end

geometry=struct('nr',nr,'nr_inner',nr_inner,'nr_outer',nr_outer,'nr_down',nr_down, ...
    'nt',nt,'nt_inner',nt_inner,'nt1',nt1,'nt2',nt2,'t_min',t_min,'t_max',t_max,'alpha',alpha,'Length',Length,...
    't_boundary',t_boundary,'tr_max',tr_max,'tl_min',tl_min);%,'r_boundary',r_boundary);


eq.R=R_new;
eq.Z=Z_new;
eq.theta=theta;
eq.rho=rho;
eq.npsi=(ones(nt,1)*rho.^2)';
eq.psi=eq.npsi*psib+psi_axis;
eq.geometry=geometry;
if icheck==1
    plot_contour(x_new,y_new,0,-1,geometry,[])
end
end





%%
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
