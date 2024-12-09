function xnew = x_adjustment(x,xend,x0)
%xnew = x_adjustment(x,xend,x0)
%   x      一维网格坐标函数
%   xend   网格的左右端点(网格的左右端点已经在x中，即x(1)和x(end)
%   x0     要求网格必须包含的点,不必写出坐标网格的端点
x0_real=zeros(size(x0));
index=zeros(size(x0));
xnew=zeros(size(x));
for i=1:length(x0)
[~,index(i)]=min(abs(x-x0(i)));
x0_real(i)=x(index(i));
end
for i=1:length(x0)-1
    xnew(index(i):index(i+1)-1)=(x(index(i):index(i+1)-1)-x0_real(i))*(x0(i+1)-x0(i))/(x0_real(i+1)-x0_real(i))+x0(i);
end
xnew(1:index(1)-1)=(x(1:index(1)-1)-xend(1))*(x0(1)-xend(1))/(x0_real(1)-xend(1))+xend(1);
xnew(index(end):end)=(x(index(end):end)-xend(2))*(x0(end)-xend(2))/(x0_real(end)-xend(2))+xend(2);
end