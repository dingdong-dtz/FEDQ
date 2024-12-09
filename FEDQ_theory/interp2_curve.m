function [xnew,ynew] = interp2_curve(x,y,method,number)
%interp2_curve(x,y,metherd,number)
%   (x,y)为一条平面连续曲线
%method表示插值方式，和interp1函数的相同
%number表示插值后曲线元素数目
s=sqrt(diff(x).^2+diff(y).^2);
[row,lie]=size(x);
if row==1
    s=[0 cumsum(s)];
else
    s=[0; cumsum(s)];
end
snew=linspace(0,s(end),number);
xnew=interp1(s,x,snew,method);
ynew=interp1(s,y,snew,method);
end