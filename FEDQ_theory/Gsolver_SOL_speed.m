function eq=Gsolver_SOL_speed(eq,nd)
%*************************************************
%  G-S euqation
%       L(psi) = Js
%  L = cJ^-1 [d_r(cJ*grr*d_r) + d_t(cJ*gtt d_t)
%         d_r(cJ*grt*d_t) + d_t(cJ*grt*d_r) ]
%    = [grr*d_r^2 + gtt*d_t^2  + 2*grt*d_r d_t]
%       +cJ^-1 {[d_r(crr)+d_t(crt)]*d_r
%              +[d_r(crt)+d_t(ctt)]*d_t}
%  Js =  - [R^2 * mu0 * P' +  FF']
%   crr=cJ*grr;
%   ctt=cJ*gtt;
%   crt=cJ*grt;
%
%*************************************************
%------------------------------ GS equation inputs
geometry=eq.geometry;
nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nt2=geometry.nt2;
nt=geometry.nt;
nr_inner=geometry.nr_inner;
nr_outer=geometry.nr_outer;
nr_down=geometry.nr_down;
nr=geometry.nr;
t_max=geometry.t_max;
t_min=geometry.t_min;
R=eq.R;  Z=eq.Z; rho=eq.rho; theta=eq.theta;
grr=eq.M.grr;        grt=eq.M.grt;    gtt=eq.M.gtt;
dcrrdr=eq.M.dcrrdr;  dcrtdt=eq.M.dcrtdt;
dcttdt=eq.M.dcttdt;  dcrtdr=eq.M.dcrtdr;
cJ=eq.M.cJ;          JupR2=eq.M.JupR2;
Js=eq.Js;
%------------------------------ finite difference method
%------------- L
L.dx2=grr;   %   grr* d_r^2
L.dy2=gtt;   %   gtt* d_t^2
L.dx1=(dcrrdr+dcrtdt).*JupR2;
L.dy1=(dcttdt+dcrtdr).*JupR2;
L.dxdy=2.0*grt;  % 2*grt*d_r d_t
%---------------
eq.L=L;
%------------ outer boundary
pre_left=sqrt((eq.boundary.Rdiv_left-eq.boundary.Rdiv_left(1)).^2+(eq.boundary.Zdiv_left-eq.boundary.Zdiv_left(1)).^2);
pre_right=sqrt((eq.boundary.Rdiv_right-eq.boundary.Rdiv_right(1)).^2+(eq.boundary.Zdiv_right-eq.boundary.Zdiv_right(1)).^2);
dgrid_left=sqrt((eq.R(nr_down:nr,t_max)-eq.R(nr,t_max)).^2+(eq.Z(nr_down:nr,t_max)-eq.Z(nr,t_max)).^2);
dgrid_right=sqrt((eq.R(nr_down:nr,t_min)-eq.R(nr,t_min)).^2+(eq.Z(nr_down:nr,t_min)-eq.Z(nr,t_min)).^2);
eq.psi(nr_down:nr,t_max)=interp1(pre_left,eq.boundary.psi_left,dgrid_left);
eq.psi(nr_down:nr,t_min)=interp1(pre_right,eq.boundary.psi_right,dgrid_right);
eq.psi(nr,t_min:t_max)=eq.boundary.psi_lcs;
eq.psi(nr_down,[t_min:nt1,nt1+nt_inner+1:t_max])=eq.boundary.psi_pri;
bdy=eq.psi(:,nt1+ceil(nt_inner/2));
%----------------- calling oprt for solving the GS equation
disp('build Matrix for GS equation ... ')
[eq.psi,eq.A,eq.B]=nabla2d_all_speed(rho,theta,L,Js,nd,bdy,geometry,R,Z,eq.psi);
disp('... GS equation solved!')
%------------------------

end