function xnew = x_adjustment(x,xend,x0)
%xnew = x_adjustment(x,xend,x0)
%X one-dimensional grid coordinate function
%The left and right endpoints of the xend grid (the left and right endpoints of the grid are already in x, i.e. x (1) and x (end))
%X0 requires that the grid must contain points, and there is no need to write the endpoints of the coordinate gridx0_real=zeros(size(x0));
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