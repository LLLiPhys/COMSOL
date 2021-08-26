function BPbands
clear all;
close all;
global Ec Ev xic xiv etac etav alpha beta
constants;
inputargs;
kk=-2:0.01:2;
[E1,H1]=bandstructure(kk,'[100]');
[E2,H2]=bandstructure(kk,'[110]');
[E3,H3]=bandstructure(kk,'[010]');

figure(1)
plot(kk,E1,kk,H1,'LineWidth',3,'Color','b');
hold on;
plot(kk,E2,kk,H2,'LineWidth',3,'Color','g');
hold on;
plot(kk,E3,kk,H3,'LineWidth',3,'Color','r');
hold off;
axis([-0.5,0.5,-2,2]);axis square;
xlabel('k [1/nm]','Fontname','Times New Roman','FontSize',36);
ylabel('E [eV]','Fontname','Times New Roman','FontSize',36);
set(gca,'LineWidth',3,'Fontname','Times New Roman','FontSize',36);

[kx,ky]=meshgrid(kk,kk);
A=Ec+xic*kx.^2+etac*ky.^2;
B=Ev-xiv*kx.^2-etav*ky.^2;
C=alpha*kx+beta*ky.^2;
Ee=(A+B+sqrt((A-B).^2+4*C.^2))/2;
Eh=(A+B-sqrt((A-B).^2+4*C.^2))/2;

figure(2)
mesh(kx,ky,Ee);
axis square;box on;grid off;

figure(3)
contour(kx,ky,Ee,10);
axis([-2,2,-2,2]);axis square;
figure(4)
surfc(kx,ky,Ee);
axis square;box on;grid off;
end

function [Ee,Eh]=bandstructure(kk,direc)
global Ec Ev xic xiv etac etav alpha beta
switch direc
    case '[100]'
        kx=kk;
        ky=0;
    case '[110]'
        kx=kk/sqrt(2);
        ky=kk/sqrt(2);
    case '[010]'
        kx=0;
        ky=kk;
end 
A=Ec+xic*kx.^2+etac*ky.^2
B=Ev-xiv*kx.^2-etav*ky.^2
C=alpha*kx+beta*ky.^2
Ee=(A+B+sqrt((A-B).^2+4*C.^2))/2;
Eh=(A+B-sqrt((A-B).^2+4*C.^2))/2;
end