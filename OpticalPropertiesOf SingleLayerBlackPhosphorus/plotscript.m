% X(N), Y(M,N): plotting data
% LW: Linewidth; FS: Fontsize
figure;
lw=4; % Linewidth
fs=40; % Fontsize
x=0:0.01:2;
y=x.*exp(-x.^2);
u=x;
v=x.^2.*exp(-x.^2);
plot(x,y,'linewidth',lw,'color','b');
hold on; 
plot(u,v,'linewidth',lw,'color','r');
hold off;
set(gca,'xlim',[0 2],'xtick',[0:0.5:2]);
set(gca,'ylim',[0 0.6],'ytick',[0:0.2:0.6]);
set(gca,'linewidth',lw,'fontname','times new roman','fontsize',fs)
xhd=xlabel('$$\hbar\omega$$ (eV)');set(xhd,'interpret','latex');
yhd=ylabel('$$\alpha(\omega)$$ (cm$$^{-1}$$)');set(yhd,'interpret','latex');
set(xhd,'fontsize',fs) %,'fontname','times new roman');
set(yhd,'fontsize',fs) %,'fontname','times new roman');
posx=get(xhd,'position');posy=get(yhd,'position');
%set(xhd,'position',[posx(1) posx(2)-0.01 posx(3)]);
%set(yhd,'position',[posy(1)-0.01 posy(2) posy(3)]);

fid1=fopen('OpticalConductivity.dat','w');
fid2=fopen('PolarizationRatio.dat','w');
fprintf(fid1,'%d %d %d %d %d\n',ET,Txx,Tyy);
fprintf(fid2,'%d %d %d\n',ET,PT);
fclose(fid1);
fclose(fid2);





