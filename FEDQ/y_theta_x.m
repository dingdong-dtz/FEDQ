function y = y_theta_x(theta,x)
%y = y_theta_x(theta,x)
%   该函数根据x和θ来确定y的值,注意函数不能在包含纵垂线的情况下使用，因为那时y不受x和θ决定
index=theta>=pi/2&theta<=pi/2*3;
index_down=theta>2*pi|theta<0;
y(~index)=-atan(x(~index).*cot(theta(~index))*pi/2)*2/pi-2*index_down(~index);
y(index)=-x(index).*cot(theta(index));
index_0=(theta==0)|(theta==2*pi);
y(index_0)=-1;
end