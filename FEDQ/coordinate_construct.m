function theta = coordinate_construct(x,y)
%theta = coordinate_construct(x,y)
%求出（x，y）对应的theta坐标值
theta=zeros(size(x));

index=find(y>-1&y<0&x>0);
[theta_x_up,~]=cart2pol(pi/2*x(index),tan(pi/2*y(index)));
theta(index)=theta_x_up+pi/2;
index=find(y>-1&y<0&x<0);
[theta_x_up,~]=cart2pol(pi/2*x(index),tan(pi/2*y(index)));
theta(index)=theta_x_up+pi/2+2*pi;
index=find(y<-1);
[theta_x_up,~]=cart2pol(pi/2*x(index),tan(pi/2*(2+y(index))));
theta(index)=theta_x_up-pi/2+(x(index)<0)*2*pi;

index=find(y>=0);
[theta_x_up,~]=cart2pol(x(index),y(index));
theta(index)=theta_x_up+pi/2;

index=x==0&y>-1&y<0;
theta(index)=0;
index=x==0&y<=-1;
theta(index)=0;
index=x>0&y==-1;
theta(index)=0;
index=x<=0&y==-1;
theta(index)=2*pi;
end