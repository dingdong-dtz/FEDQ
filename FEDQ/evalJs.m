function Js=evalJs(gpsi,Pprime,FFprime,gpsi_pri,Pprime_pri,FFprime_pri,R,psi,geometry)
%*****************************************************
% function Js=evalJs(gpsi,Pprime,FFprime,R,nrho,npsi)
%  evaluate right hand side of GS equaiton
%       Js = - [ R^2 * mu0 * P' + FF'];
% inputs:
%       gpsi,       normalized psi_p grids for P' and FF'
%       Pprime,     P'
%       FFprime     FF'
%       gpsi_pri
%       Pprime_pri
%       FFprime_pri
%R, major radius of the surfaces (arranged in matrix position)
%Psi, psi_p (arranged by matrix position) psi should correspond to gpsi, either normalized or unnormalized
%******************************************************
nt1=geometry.nt1;
nt_inner=geometry.nt_inner;
nt2=geometry.nt2;
nt=geometry.nt;
nr_inner=geometry.nr_inner;
nr_outer=geometry.nr_outer;
nr_down=geometry.nr_down;
nr=geometry.nr;
t_min=geometry.t_min;
t_max=geometry.t_max;
Pprime_new=zeros(nr,nt);
FFprime_new=zeros(nr,nt);
Pprime_new(1:nr_inner,nt1+1:nt1+nt_inner)=interp1(gpsi,Pprime,psi(1:nr_inner,nt1+1:nt1+nt_inner),'pchip');  % Pprime is defined in npsi
FFprime_new(1:nr_inner,nt1+1:nt1+nt_inner)=interp1(gpsi,FFprime,psi(1:nr_inner,nt1+1:nt1+nt_inner),'pchip');
Pprime_new(nr_inner+1:end,:)=interp1(gpsi,Pprime,psi(nr_inner+1:end,:),'pchip');  % Pprime is defined in npsi
FFprime_new(nr_inner+1:end,:)=interp1(gpsi,FFprime,psi(nr_inner+1:end,:),'pchip');
Pprime_new(1:nr_inner,[1:nt1 nt1+nt_inner+1:nt])=interp1(gpsi_pri,Pprime_pri,psi(1:nr_inner,[1:nt1 nt1+nt_inner+1:nt]),'pchip');  % Pprime is defined in npsi
FFprime_new(1:nr_inner,[1:nt1 nt1+nt_inner+1:nt])=interp1(gpsi_pri,FFprime_pri,psi(1:nr_inner,[1:nt1 nt1+nt_inner+1:nt]),'pchip');
Js = -[R.^2.*Pprime_new + FFprime_new];
end
