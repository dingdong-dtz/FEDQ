function [x,y] = rotate(R,Z,R0,Z0,alpha,L)
%[x,y] = rotate(R,Z,R0,Z0,alpha,L)
%将坐标系顺时针旋转pi/2-α,以（R0,Z0）为中心并缩小L倍
x=((R-R0)*sin(alpha)-(Z-Z0)*cos(alpha))/L;
y=((R-R0)*cos(alpha)+(Z-Z0)*sin(alpha))/L;
end