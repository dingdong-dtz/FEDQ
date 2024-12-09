function [R0,Z0,dpsi_dr,psi_axis_fit]=fit_axis_initial(R,Z,psi,R0_pre,Z0_pre,axis_range,icheck)
%**************************************
% function [R0,Z0,npsi,psibm,psi_axis,sIp,r0,t0]=findaxis(eq,icheck)
% find magnetic axis
% psi_axis is the smaller one within 'psi_axis_fit' and the smallest number of all psi
% 磁轴处等势面接近椭圆形，用内层磁面拟合确定磁轴位置
%*************************************
index=((R-R0_pre).^2+(Z-Z0_pre).^2)<axis_range^2;
near_R=R(index);
near_Z=Z(index);
near_psi=psi(index);

roundEqn = '(a*(x-R0).^2+b*(y-Z0).^2)+psi_axis';
startPoints = [R0_pre,Z0_pre,1,1,0];%利用已知条件确定初始值
[f1,gof] = fit([near_R(:),near_Z(:)],near_psi(:),roundEqn,'Start', startPoints);
%coeffnames(f1)
R0=f1.R0;
Z0=f1.Z0;
dpsi_dr=[f1.a f1.b];
psi_axis_fit=f1.psi_axis;
if icheck
    disp('均方根误差（标准误差）')
    gof.rmse
    figure
    hold on
    scatter(R0,Z0);
end
if ((R0-R0_pre)^2+(Z0-Z0_pre)^2)>axis_range^2/9
[R0,Z0,dpsi_dr,psi_axis_fit]=fit_axis_initial(R,Z,psi,R0,Z0,axis_range,icheck);
end
end