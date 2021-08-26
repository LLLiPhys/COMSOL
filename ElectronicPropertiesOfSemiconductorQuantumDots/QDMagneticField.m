%%% Main program
% Model confined electron and hole energy states 
% in a parabolic quantum dot in a magnetic field 
% by solving the effective-mass equation with 
% appropriate boundary conditions using COMSOL

clear all;
close all;
% Fundamental physical constants
hbar=1.054571628e-34; % [J*s]
m0=9.10938215e-31; % [Kg]
q=1.602176487e-19; % [C]
% Scaling parameters in the Schrodinger equation
alfa0=hbar^2/m0/1^2/1e-18/q*1e3;
hw0=hbar^2/m0/1^2/1e-18/q*1e3;
hwc=hbar*q*1/m0/q*1e3;
% Boundary condiction type
bctype='DBC'; % Dirichlet type (others: Neumann & Periodic)
% Mesh grid points along the radial direction
nx=101;x=linspace(0,1,nx);dx=x(2)-x(1);
% Step length for varying the magnetic field
dB=0.5;Bmin=0;Bmax=20;
% Angular and radial quantum numbers
mmax=5;nmax=5;
% Sweeping parameters for the model calculation
Blist=0:dB:Bmax;
mlist=-mmax:mmax;
nlist=0:nmax-1;
nB=length(Blist);
nm=length(mlist);
nn=length(nlist);
% Material parameters for the quantum disk
R=10; % Radius of the quantum disk [nm]
mes=0.067; % Electron effective mass [m0]

%hw0e=hbar^2/me/R^2/1e-18/q*1e3;
%hwce=hbar*q*B/me/q*1e3;
%alfae=hbar^2/me/R^2/1e-18/q*1e3;
%betae=hwce/alfae;
%gamae=hw0e/alfae;

me=zeros(nB,nm,nn);ne=zeros(nB,nm,nn);ee=zeros(nB,nm,nn);ue=zeros(nB,nm,nn,nx);
for iB=1:nB
    B=Blist(iB);
    hw0e=hw0/mes/R^2;
    hwce=hwc*B/mes;
    alfae=alfa0/mes/R^2;
    betae=hwce/alfae;
    gamae=hw0e/alfae;
    for im=1:nm
        m=mlist(im);
        feme=[];lbdae=[];ene=[];elist=[];eindex=[];
        %femh=[];lbdah=[];enh=[];hlist=[];hindex=[];
        [feme]=femeigsolver(alfae,betae,gamae,m);
        lbdae=feme.sol.lambda;lbdae=real(lbdae);
        ene=sort(lbdae,'ascend');elist=ene(find(ene>0,nmax,'first'));
        %enh=sort(lbdah,'descend');hlist=enh(find(enh<0,nmax,'first'));
        eindex=[];
        %hindex=[];
        for ie=1:length(elist)
            eindex=[eindex find(lbdae==elist(ie))];
        end
        %for ih=1:length(hlist)
            %hindex=[hindex find(lbdah==hlist(ih))];
        %end
        elist=[];
        %hlist=[];
        elist=lbdae(eindex);
        %hlist=lbdah(hindex);
        fprintf('===============B=%3.1f/m=%0.0f===============\n',B,m); 
        disp('=====Electron Energy (eV)=======Hole Energy (eV)=====');
        disp([elist']);
        intue=[];
        %intuh=[];
        intue=postint(feme,'u*conj(u)*x','solnum',eindex);
        %intuh=postint(femh,'u*conj(u)*x','solnum',hindex);
        fe=[];
        %fh=[];
        fe=postinterp(feme,'u',x,'solnum',eindex);
        %fh=postinterp(femh,'u',x,'solnum',hindex);
        for ie=1:length(intue)
            fe(ie,:)=fe(ie,:)./sqrt(intue(ie))./R;
        end
        %for ih=1:length(intuh)
            %fh(ih,:)=fh(ih,:)./sqrt(intuh(ih))./R;
        %end
        for in=1:nn
            n=nlist(in);
            me(iB,im,in)=m;
            %mh(iB,im,in)=m;
            ne(iB,im,in)=n;
            %nh(iB,im,in)=n;
            ee(iB,im,in)=elist(in);
            %eh(iB,im,in)=hlist(in);
            ue(iB,im,in,1:nx)=ue(in,1:nx);
            %uh(iB,im,in,1:nx)=uh(in,1:nx);
        end
        clear intue fe 
        %clear intuh fh
        clear feme lbdae ene elist eindex
        %clear femh lbdah enh hlist hindex
    end
end

Ee=zeros(nB,nm,nn);
for iB=1:nB
    B=Blist(iB);
    hw0e=hw0/mes/R^2;
    hwce=hwc*B/mes;
    hwe=sqrt(hwce^2+4*hw0e^2);
    for im=1:nm
        m=mlist(im);
        for in=1:nn
            n=nlist(in);
            Ee(iB,im,in)=hwe*(n+(abs(m)+1)/2)+hwce*m/2;
        end
    end
end

B=zeros(nB,1);B(1:nB)=Blist(1:nB);
Emin=0;Emax=50;


%%% Eigensolver using FEM in COMSOL

function [fem]=femeigsolver(alfa,beta,gama,m)

flclear fem

% Constants
fem.const = {'alfa',alfa,'beta',beta,'gama',gama,'m',m};

% Geometry
g1=solid1([0,10]);

% Analyzed geometry
clear s
s.objs={g1};
s.name={'I1'};
s.tags={'g1'};

fem.draw=struct('s',s);
fem.geom=geomcsg(fem);

% Initialize mesh
fem.mesh=meshinit(fem,'hmax',0.01);

% (Default values are not included)

% Application mode 1
clear appl
appl.mode.class = 'FlPDEC';
appl.sshape = 2;
appl.assignsuffix = '_c';
clear bnd

bnd.type = 'dir';
bnd.h = {1,0};
bnd.ind = [2,1];
appl.bnd = bnd;

% bnd.type = {'dir','neu'};
% bnd.h = {0,1};
% bnd.ind = [1,2];
% appl.bnd = bnd;

clear equ
equ.f = 0;
equ.c = 'alfa/2';
equ.a = 'alfa*m^2/2/x^2+alfa*beta*m/2+alfa*(gama^2+beta^2/4)*x^2/2';
equ.be = '-alfa/2/x';
equ.ind = [1];
appl.equ = equ;
fem.appl{1} = appl;
fem.frame = {'ref'};
fem.border = 1;

% ODE Settings
clear ode
clear units;
units.basesystem = 'SI';
ode.units = units;
fem.ode=ode;

% Multiphysics
fem=multiphysics(fem);

% Extend mesh
fem.xmesh=meshextend(fem);

% Solve problem
fem.sol=femeig(fem,'solcomp',{'u'},'outcomp',{'u'},'blocksize','auto','neigs',50);

% Save current fem structure for restart purposes
%fem0=fem;

end
plot(B,ee(:,:),'r-',B,Ee(:,:),'b--','linewidth',2);
axis square;axis([Bmin Bmax Emin Emax]);

