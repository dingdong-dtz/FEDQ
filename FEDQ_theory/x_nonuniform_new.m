function xn=x_nonuniform_new(xends,nx,x0,d,rate,icheck)
%**********************************************************
%  xn=x_nonuniform(xends,nx,x0,d,rate,icheck)
%  generate smooth non-uniform grids xn 
%      with a grid number nx, 
%      in the range [xends(1), xends(2)],
%      and increase the grid density roughly 
%      in the range [x0-d x0+d] with a factor 'rate'.
%   Here x0(:), d(:), rate(:) can be arrays, 
%      and they should have the same length.
%   If icheck==1, plot the distribution of xn
%  e.g. 
%   increase grid density near one point at x0=0.8:
%   xn=x_nonuniform([0  1],65,0.8,0.2,20,0);
%   increase grid density near two points at x0=[0.5,1]:
%   xn=x_nonuniform([0  1],65,[0.5 1],[0.01 0.03],[10 20],0);
%   improve by dtz
%**********************************************************
%----- initial uniform profiles
%----- continous transform function
%--generate a function with steep gradient near x0
x=linspace(xends(1),xends(2),nx);
xpre=linspace(xends(1),xends(2),nx*100);
f=ones(1,nx*100);
for i=1:length(x0)
    f=f+rate(i)*exp(-(xpre-x0(i)).^2/d(i)^2);
end
f=cumsum(f);
for i=1:length(x0)
    f=f-0.5*rate(i)*exp(-(xpre-x0(i)).^2/d(i)^2);
end
%-- normalization to range of x
f=f*(xends(2)-xends(1))/(f(end)-f(1));
f=f-f(1)+xends(1);
f(end)=xends(2);
%----- new profiles with non-uniform grids
xn=interp1(f,xpre,x,'pchip');
xn(end)=xends(2);
xn(1)  =xends(1);
%----- plot results
if icheck==1
    figure('Name','meshgrid_choose')
    clf
    plot(xn,xn*0,'*')
end
end