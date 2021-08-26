function Ef=fermilevel(T,N,A,B)
global Ec Ev xic xiv etac etav alpha beta;
global kT Ns;
inputargs;
kT=0.026*T/300;
Ns=N*1e-16;
Ef=fzero(@rtsfermi,[A B]);
end
function f=rtsfermi(x)
global Ec Ev xic xiv etac etav alpha beta;
global kT Ns;
global mu;
mu=x;
xa=0;xb=inf;ya=0;yb=2*pi;
Q=integral2(@intfermi,xa,xb,ya,yb);
f=Q/2/pi^2/Ns-1.0;
end
function f=intfermi(x,y)
global Ec Ev xic xiv etac etav alpha beta;
global kT Ns;
global mu;
u=cos(y);
v=sin(y);
A=Ec+xic*x.^2.*u.^2+etac*x.^2.*v.^2;
B=Ev-xiv*x.^2.*u.^2-etav*x.^2.*v.^2;
C=alpha*x.*u+beta*x.^2.*v.^2;
Exy=(A+B+sqrt((A-B).^2+4*C.^2))/2;
z=(Exy-mu)./kT;
f=1./(exp(z)+1).*x;
%f=feval(@fermi,z).*x;
end
function f=fermi(x)
f=1.0.*(x<=-200)+1./(exp(x)+1).*(abs(x)<200)+0.0.*(x>=200);
end