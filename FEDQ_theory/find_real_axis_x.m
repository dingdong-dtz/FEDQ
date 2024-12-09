function [R0_real,Z0_real,psi0_theory,Rx_real,Zx_real,psi_X_theory] = find_real_axis_x(R0,Z0,Rx,Zx)
%[R0_real,Z0_real,psi0_theory,Rx_real,Zx_real,psi_X_theory] = find_real_axis_x(R0,Z0,Rx,Zx)
%This function calculates the magnetic axis position and x-point position for a given analytical example, which is used to calculate the error in program calculations
%R0, Z0, Rx, Zx are estimated positions of the speculated magnetic axis and x-point
options = optimoptions(@fminunc,'OptimalityTolerance',1e-10);
fun=@(x) theory_sample(x(1),x(2));
[axis_theory,psi0_theory] = fminunc(fun,[R0,Z0],options);
R0_real=axis_theory(1); Z0_real=axis_theory(2);
%Find the position of point X by taking the zero point of the derivative function
syms RR  ZZ
psi=theory_sample(RR,ZZ);
psi_R=diff(psi,RR);
psi_Z=diff(psi,ZZ);
%gpp_theory=matlabFunction(psi_R^2+psi_Z^2);
dpsidR=matlabFunction(psi_R);
dpsidZ=matlabFunction(psi_Z);
grad=@(x) [dpsidR(x(1),x(2));dpsidZ(x(1),x(2))];
options = optimoptions(@fsolve,'OptimalityTolerance',1e-30,'Display','iter');
[X_theory,~] = fsolve(grad,[Rx,Zx],options);
Rx_real=X_theory(1);Zx_real=X_theory(2);
psi_X_theory=theory_sample(Rx_real,Zx_real);
end