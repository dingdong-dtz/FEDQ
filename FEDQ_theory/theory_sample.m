function psi = theory_sample(R,Z)
R1=5;
R2=7;
Z2=4.5;
A=-0.0005;
B=0.01;
psi=(1-R.^2/R2^2-Z.^2/Z2^2).*(R.^2-R1^2)+A*(3*R.^2-4*Z.^2).*R.^2.*Z+B*R.^2;
psi=-psi;
end