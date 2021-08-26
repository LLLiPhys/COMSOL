function BPoptics
clear all;
close all;
global Ec Ev xic xiv etac etav alpha beta;
inputargs;
T=77;N=1e12;del=1e-3;
Ef=fermilevel(T,N,-5,5);
fprintf(1,'T= %d N = %d Ef = %d \n',T,N,Ef);
% Interband optical conductivity tensor:
% J_xx,J_yy, and J_xy
% Real and imaginary parts of interband sigma
% Intraband optical conductivity tensor:
% L_xx,L_yy, and L_xy
% Real and imaginary parts of intraband sigma
EJ=[];EL=[];ET=[];PJ=[];PL=[];PT=[];
J1=[];J2=[];J3=[];L1=[];L2=[];L3=[];T1=[];T2=[];T3=[];
DelJ=[];DelL=[];DelT=[];
Jxx=[];Jyy=[];Jxy=[];Lxx=[];Lyy=[];Lxy=[];Txx=[];Tyy=[];Txy=[];
fid1=fopen('5KOC.dat','w');gid1=fopen('5KPR.dat','w');
fid2=fopen('77KOC.dat','w');gid2=fopen('77KPR.dat','w');
fid3=fopen('150KOC.dat','w');gid3=fopen('150KPR.dat','w');
fid4=fopen('300KOC.dat','w');gid4=fopen('300KPR.dat','w');
for IT=1:4
    if IT==1;T=5;end
    if IT==2;T=77;end
    if IT==3;T=150;end
    if IT==4;T=300;end
    Ef=fermilevel(T,N,-5,5);
    fprintf(1,'T= %d N = %d Ef = %d \n',T,N,Ef);
for hw=0:1e-3:3 
    [J1,J2,J3]=intersigma(T,Ef,del,hw);
    [L1,L2,L3]=intrasigma(T,Ef,del,hw);
    if J1<1e-3;J1=1e-3;end
    if J2<1e-3;J2=1e-3;end
    if L1<1e-3;L1=1e-3;end
    if L2<1e-3;L2=1e-3;end
    T1=J1+L1;T2=J2+L2;T3=J3+L3;
    DelJ=abs((J2-J1)./(J2+J1));
    DelL=abs((L2-L1)./(L2+L1));
    DelT=abs((T2-T1)./(T2+T1));
    EJ=[EJ;hw];Jxx=[Jxx;J1];Jyy=[Jyy;J2];Jxy=[Jxy;J3];
    EL=[EL;hw];Lxx=[Lxx;L1];Lyy=[Lyy;L2];Lxy=[Lxy;L3];
    ET=[ET;hw];Txx=[Txx;T1];Tyy=[Tyy;T2];Txy=[Txy;T3];
    PJ=[PJ;DelJ];PL=[PL;DelL];PT=[PT;DelT];
    fprintf(1,'hw= %d ReJxx = %d ReJyy = %d \n',hw,J1,J2);
    fprintf(1,'hw= %d ReLxx = %d ReLyy = %d \n',hw,L1,L2);
    if T==5;
    fprintf(fid1,'%d %d %d %d %d %d %d\n',hw,T1,T2,L1,L2,J1,J2);
    fprintf(gid1,'%d %d %d %d\n',hw,DelT,DelL,DelJ);
    end
    if T==77;
    fprintf(fid2,'%d %d %d %d %d %d %d\n',hw,T1,T2,L1,L2,J1,J2);
    fprintf(gid2,'%d %d %d %d\n',hw,DelT,DelL,DelJ);
    end
    if T==150;
    fprintf(fid3,'%d %d %d %d %d %d %d\n',hw,T1,T2,L1,L2,J1,J2);
    fprintf(gid3,'%d %d %d %d\n',hw,DelT,DelL,DelJ);
    end
    if T==300;
    fprintf(fid4,'%d %d %d %d %d %d %d\n',hw,T1,T2,L1,L2,J1,J2);
    fprintf(gid4,'%d %d %d %d\n',hw,DelT,DelL,DelJ);
    end
end
end
fclose(fid1);fclose(gid1);
fclose(fid2);fclose(gid2);
fclose(fid3);fclose(gid3);
fclose(fid4);fclose(gid4);

subplot(1,2,1);
h1=plot(ET,Jxx,'b',ET,Jxy,'r');set(h1,'linewidth',3);
set(gca,'linewidth',3,'fontname','times new roman','fontsize',36);
xh=xlabel('$$\hbar\omega$$ (eV)');set(xh,'interpret','latex');
yh=ylabel('$$Re[\sigma/\sigma_0]$$');set(yh,'interpret','latex');
set(xh,'fontsize',36);set(yh,'fontsize',36) %('fontname','times new roman');
subplot(1,2,2);
h2=plot(ET,PT);set(h2,'linewidth',3);
set(gca,'linewidth',3,'fontname','times new roman','fontsize',36);
xh=xlabel('$$\hbar\omega$$ (eV)');set(xh,'interpret','latex');
yh=ylabel('$$Polarization Ratio$$');set(yh,'interpret','latex');
set(xh,'fontsize',36);set(yh,'fontsize',36) %('fontname','times new roman');
end
function Ef=fermilevel(T,N,A,B)
global Ec Ev xic xiv etac etav alpha beta;
global kT Ns;
kT=0.026*T/300;Ns=N*1e-16;
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
Ee=(A+B+sqrt((A-B).^2+4*C.^2))/2;
ze=(Ee-mu)/kT;
%f=1./(exp(ze)+1).*x;
f=feval(@fermi,ze).*x;
end

function [J1,J2,J3]=intersigma(T0,Ef0,del0,hw0)
global Ec Ev xic xiv etac etav alpha beta;
global T Ef del hw;
global flag
T=T0;Ef=Ef0;del=del0;hw=hw0;
xa=0;xb=inf;ya=0;yb=2*pi;
J1=0;J2=0;J3=0;
flag='ReJxx';P=integral2(@intinterband,xa,xb,ya,yb);J1=P/2/pi^2;
% flag='ImJxx';P=integral2(@intinterband,xa,xb,ya,yb);J1(2)=P/2/pi^2;
flag='ReJyy';P=integral2(@intinterband,xa,xb,ya,yb);J2=P/2/pi^2;
% flag='ImJyy';P=integral2(@intinterband,xa,xb,ya,yb);J2(2)=P/2/pi^2;
flag='ReJxy';P=integral2(@intinterband,xa,xb,ya,yb);J3=P/2/pi^2;
%flag='ImJxy';P=integral2(@intinterband,xa,xb,ya,yb);J3(2)=P/2/pi^2;
end

function [L1,L2,L3]=intrasigma(T0,Ef0,del0,hw0)
global Ec Ev xic xiv etac etav alpha beta;
global T Ef del hw;
global flag
T=T0;Ef=Ef0;del=del0;hw=hw0;
xa=0;xb=inf;ya=0;yb=2*pi;
L1=0;L2=0;L3=0;
flag='ReKxx';P=integral2(@intintraband,xa,xb,ya,yb);L1=P/2/pi^2;
%flag='ImKxx';P=integral2(@intintraband,xa,xb,ya,yb);L1(2)=P/2/pi^2;
flag='ReKyy';P=integral2(@intintraband,xa,xb,ya,yb);L2=P/2/pi^2;
% flag='ImKyy';P=integral2(@intintraband,xa,xb,ya,yb);L2(2)=P/2/pi^2;
% flag='ReKxy';P=integral2(@intintraband,xa,xb,ya,yb);L3(1)=P/2/pi^2;
% flag='ImKxy';P=integral2(@intintraband,xa,xb,ya,yb);L3(2)=P/2/pi^2;
end

function f=intinterband(x,y)
global Ec Ev xic xiv etac etav alpha beta;
global T Ef del hw;
global flag;
kT=0.026*T/300;
u=cos(y);
v=sin(y);
A=Ec+xic*x.^2.*u.^2+etac*x.^2.*v.^2;
B=Ev-xiv*x.^2.*u.^2-etav*x.^2.*v.^2;
C=alpha*x.*u+beta*x.^2.*v.^2;
D=((A-B)+sqrt((A-B).^2+4*C.^2))/2;
Ax=2*xic*x.*u;Ay=2*etac*x.*v;
Bx=-2*xiv*x.*u;By=-2*etav*x.*v;
Cx=alpha;Cy=2*beta*x.*v;
Ee=(A+B+sqrt((A-B).^2+4*C.^2))/2;
Eh=(A+B-sqrt((A-B).^2+4*C.^2))/2;
Gxx=(C.*D.*(Ax-Bx)+(C-D).*(C+D).*Cx).^2./(C.^2+D.^2).^2;
Gyy=(C.*D.*(Ay-By)+(C-D).*(C+D).*Cy).^2./(C.^2+D.^2).^2;
Gxy=(C.*D.*(Ax-Bx)+(C-D).*(C+D).*Cx).*(C.*D.*(Ay-By)+(C-D).*(C+D).*Cy)./(C.^2+D.^2).^2;
%ze=(Ee-Ef)./kT;zh=(Eh-Ef)./kT;fe=1./(exp(ze)+1);fh=1./(exp(zh)+1);
ze=(Ee-Ef)./kT;zh=(Eh-Ef)./kT;fe=feval(@fermi,ze);fh=feval(@fermi,zh);
switch flag
    case 'ReJxx'
        f=(fh-fe)./(Ee-Eh).*Gxx.*del./((Eh-Ee+hw).^2+del.^2).*x;
    case 'ImJxx'
        f=(fh-fe)./(Ee-Eh).*Gxx.*hw./((Eh-Ee+hw).^2+del.^2).*x;
     case 'ReJyy'
        f=(fh-fe)./(Ee-Eh).*Gyy.*del./((Eh-Ee+hw).^2+del.^2).*x;
    case 'ImJyy'
        f=(fh-fe)./(Ee-Eh).*Gyy.*hw./((Eh-Ee+hw).^2+del.^2).*x;
    case 'ReJxy'
        f=(fh-fe)./(Ee-Eh).*Gxy.*del./((Eh-Ee+hw).^2+del.^2).*x;
    case 'ImJxy'
        f=(fh-fe)./(Ee-Eh).*Gxy.*(Eh-Ee+hw)./((Eh-Ee+hw).^2+del.^2).*x;
end     
end

function f=intintraband(x,y)
global Ec Ev xic xiv etac etav alpha beta;
global T Ef del hw;
global flag;
kT=0.026*T/300;
u=cos(y);
v=sin(y);
A=Ec+xic*x.^2.*u.^2+etac*x.^2.*v.^2;
B=Ev-xiv*x.^2.*u.^2-etav*x.^2.*v.^2;
C=alpha*x.*u+beta*x.^2.*v.^2;
D=((A-B)+sqrt((A-B).^2+4*C.^2))/2;
Ax=2*xic*x.*u;Ay=2*etac*x.*v;
Bx=-2*xiv*x.*u;By=-2*etav*x.*v;
Cx=alpha;Cy=2*beta*x.*v;
Ee=(A+B+sqrt((A-B).^2+4*C.^2))/2;
Eh=(A+B-sqrt((A-B).^2+4*C.^2))/2;
Gxx=(D.^2.*Ax+C.^2.*Bx+2*C.*D.*Cx).^2./(C.^2+D.^2).^2;
Gyy=(D.^2.*Ay+C.^2.*By+2*C.*D.*Cy).^2./(C.^2+D.^2).^2;
Gxy=(D.^2.*Ax+C.^2.*Bx+2*C.*D.*Cx).*(D.^2.*Ay+C.^2.*By+2*C.*D.*Cy)./(C.^2+D.^2).^2;
%ze=(Ee-Ef)./kT;zh=(Eh-Ef)./kT;fe=1./(exp(ze)+1);fh=1./(exp(zh)+1);
ze=(Ee-Ef)./kT;zh=(Eh-Ef)./kT;fe=feval(@fermi,ze);fh=feval(@fermi,zh);
switch flag
    case 'ReKxx'
        f=fe.*(1-fe)./kT.*Gxx.*del./(hw.^2+del.^2).*x;
    case 'ImKxx'
        f=fe.*(1-fe)./kT.*Gxx.*hw./(hw.^2+del.^2).*x;
     case 'ReKyy'
        f=fe.*(1-fe)./kT.*Gyy.*del./(hw.^2+del.^2).*x;
    case 'ImKyy'
        f=fe.*(1-fe)./kT.*Gyy.*hw./(hw.^2+del.^2).*x;
    case 'ReKxy'
        f=fe.*(1-fe)./kT.*Gxy.*del./(hw.^2+del.^2).*x;
    case 'ImKxy'
        f=fe.*(1-fe)./kT.*Gxy.*hw./(hw.^2+del.^2).*x;
end     
end

function f=fermi(x)
f=1.0.*(x<=-200)+1./(exp(x)+1).*(abs(x)<200)+0.0.*(x>=200);
end