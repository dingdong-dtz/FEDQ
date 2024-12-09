function [dpsidchi,dpsids]=cald_2d(psi,chi,s)
%******************************
% function [dpsidchi,dpsids]=cald_2d(psi,chi,s)
%     2D function psi(ns,nchi)
%            psi(s,chi)
%          or  psi(r,z)
%  calculate dpsi/dchi and dpsi/ds
%   for non-uniform grids
%
%     by Y. Sun
%******************************
ns=length(s); nchi=length(chi);
s=reshape(s,[1,ns]);
chi=reshape(chi,[1,nchi]);
%-----
[ls,lchi]=size(psi);
dchiu=[diff(chi) chi(lchi)-chi(lchi-1)];
dchid=[chi(2)-chi(1) diff(chi)];
dsu=[diff(s) s(ls)-s(ls-1)];
dsd=[s(2)-s(1) diff(s)];

[Dchiu,Dsu]=meshgrid(dchiu,dsu);
[Dchid,Dsd]=meshgrid(dchid,dsd);

nddpu=[psi(:,2:lchi)-psi(:,1:(lchi-1)),psi(:,lchi)-psi(:,lchi-1)];
nddpd=[psi(:,2)-psi(:,1), psi(:,2:lchi)-psi(:,1:(lchi-1))];
dpsidchi=(Dchid.*(nddpu./Dchiu)+Dchiu.*(nddpd./Dchid))./(Dchid+Dchiu);

%---second order accuracy for boundary points
[a,b,c,d]=clb(chi);
dpsidrbl=[a*psi(:,1)+b*psi(:,2)+c*psi(:,3)]/d;
[a,b,c,d]=crb(chi);
dpsidrbr=[a*psi(:,end)+b*psi(:,end-1)+c*psi(:,end-2)]/d;
dpsidchi(:,1)=dpsidrbl;
dpsidchi(:,end)=dpsidrbr;
%---

nddpu=[psi(2:ls,:)-psi(1:(ls-1),:);psi(ls,:)-psi(ls-1,:)];
nddpd=[psi(2,:)-psi(1,:); psi(2:ls,:)-psi(1:(ls-1),:)];

dpsids=(Dsd.*(nddpu./Dsu)+Dsu.*(nddpd./Dsd))./(Dsd+Dsu);

%---second order accuracy for boundary points
[a,b,c,d]=clb(s);
dpsidrbl=[a*psi(1,:)+b*psi(2,:)+c*psi(3,:)]/d;
[a,b,c,d]=crb(s);
dpsidrbr=[a*psi(end,:)+b*psi(end-1,:)+c*psi(end-2,:)]/d;
dpsids(1,:)=dpsidrbl;
dpsids(end,:)=dpsidrbr;
%---
end
function [a,b,c,d]=clb(r)
%********************************************
% coefficients for 1st order derivatives
% with second order accuracy
%    for left boundary point
%  Du= [a*U+b*U(r+h1)+c*U(r+h1+h2)]/d
%********************************************
h1=r(2)-r(1); h2=r(3)-r(2);
a=-(2*h1+h2)*h2;    b=(h1+h2)^2;  c=-h1^2;   d=h1*h2*(h1+h2);
end
function [a,b,c,d]=crb(r)
%********************************************
% coefficients for 1st order derivatives
% with second order accuracy
%    for right boundary point
%  Du= [a*U+b*U(r-h1)+c*U(r-h1-h2)]/d
%********************************************
h1=r(end)-r(end-1); h2=r(end-1)-r(end-2);
a=(2*h1+h2)*h2;    b=-(h1+h2)^2;  c=h1^2;   d=h1*h2*(h1+h2);
%---
end