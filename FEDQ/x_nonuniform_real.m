function xn=x_nonuniform_real(xends,nx,x_interval,f_interval,icheck)
%**********************************************************
%  xn=x_nonuniform_real(xends,nx,x_interval,f_interval,icheck)
%  generate non-uniform grids xn 
%      with a grid number nx, 
%      in the range [xends(1), xends(2)],
%   If icheck==1, plot the distribution of xn
%   f_interval(x_interval)是得到的网格点所要满足的间隔比例函数，以行向量输入
%**********************************************************
%----- initial uniform profiles
%----- continous transform function
%--generate a function with steep gradient near x0
f_interval=f_interval(:)';%把输入的化成行向量
f_interval=f_interval/trapz(x_interval,f_interval)*(xends(2)-xends(1));
%f_interval=max(f_interval,ones(size(f_interval))/nx/100);
x_interval=x_interval(:)';%把输入的化成行向量
L=length(x_interval);
xn=interp1(0:L-1,x_interval,linspace(0,L-1,nx));
xn([1,nx])=x_interval([1,L]);%避免边界处插值出现问题
f=interp1(x_interval,f_interval,xn);
error=1;
k=1;
while error>0.00001&&k<100
    diff_rho=cald_1d(xn,1:nx);
    fcum=cumtrapz(xn,f./diff_rho);
    fcum=fcum/fcum(end)*(xends(2)-xends(1))+xends(1);
    xn_new=interp1(xn,fcum,xn);
    xn_new([1,nx])=xn([1,nx]);%避免边界处插值出现问题
    error=max(abs((xn_new-xn)/(xends(2)-xends(1))));
    f=interp1(x_interval,f_interval,xn_new);
    % if mod(k,10)==0
    %  figure
    %  hold on
    %  plot(xn,diff_rho*f(5)/diff_rho(5),'r*-');
    %  plot(xn_new,f,'b-')
    % end
    xn=xn_new;
    k=k+1;
end
xn(end)=xends(2);
xn(1)  =xends(1);
%----- plot results
if icheck==1
    figure
    hold on
    plot(x_interval,f_interval,'b.-');
    % plot(xn,f,'b.-');
    plot(xn,diff_rho*f(5)/diff_rho(5),'r*-');
    figure('Name','meshgrid_choose')
    clf
    plot(xn,xn*0,'*')
end

% if x_interval(1)>xends(1)

% figure
% plot(x_interval_new,f_interval_new)
% else
% end
end

function [dpsidr]=cald_1d(psi,r)
%******************************
% function [dpsidr]=cald_1d(psi,r)
%  calculate dpsi/dr 
%   for non-uniform grids
%******************************
     lr=length(psi); 
     %--
     s0=size(r);
     psi=reshape(psi,[1,lr]);r=reshape(r,[1,lr]);  %  avoid errors with one collomn
     %--
     dru=[diff(r) r(lr)-r(lr-1)];
     drd=[r(2)-r(1) diff(r)];   
     nddpu=[psi(2:lr)-psi(1:(lr-1)),psi(lr)-psi(lr-1)];
     nddpd=[psi(2)-psi(1), psi(2:lr)-psi(1:(lr-1))];
     dpsidr=(drd.*(nddpu./dru)+dru.*(nddpd./drd))./(drd+dru);
    %---second order accuracy for boundary points
        [a,b,c,d]=clb(r);
        dpsidrbl=[a*psi(1)+b*psi(2)+c*psi(3)]/d;
         [a,b,c,d]=crb(r);
        dpsidrbr=[a*psi(end)+b*psi(end-1)+c*psi(end-2)]/d;       
     dpsidr(1)=dpsidrbl;
     dpsidr(end)=dpsidrbr;
     %--
     dpsidr=reshape(dpsidr,s0); %  change its shape back     
end
function [a,b,c,d]=clb(r)
%********************************************
% coefficients for 1st order derivatives
% with second order accuracy 
%    for left boundary point
%  Du= [a*U+b*U(r+h1)+c*U(r+h1+h2)]/d
%********************************************
        h1=r(2)-r(1); h2=r(3)-r(2);
        a=-(2*h1+h2)*h2;    b=(h1+h2)^2;  c=-h1^2;   d=h1*h2*(h1+h2);
end      
function [a,b,c,d]=crb(r)
%********************************************
% coefficients for 1st order derivatives
% with second order accuracy 
%    for right boundary point
%  Du= [a*U+b*U(r-h1)+c*U(r-h1-h2)]/d
%********************************************
        h1=r(end)-r(end-1); h2=r(end-1)-r(end-2);
        a=(2*h1+h2)*h2;    b=-(h1+h2)^2;  c=h1^2;   d=h1*h2*(h1+h2);
end