function x = x_theta_y(theta,y)
%x = x_theta_y(theta,y)
%   该函数根据y和θ来确定x的值
x=zeros(size(y));
index=theta>=pi/2&theta<=pi/2*3;
index_down=theta>2*pi|theta<0;
x(~index)=-2/pi*tan(theta(~index)).*tan(pi/2*(y(~index)+2*index_down(~index)));
x(index)=-y(index).*tan(theta(index));
%y(~index)=-atan(x(~index).*cot(theta(~index))*pi/2)*2/pi-2*index_down(~index);
%y(index)=-x(index).*cot(theta(index));
%index_0=theta==0;
%y(index_0)=
end