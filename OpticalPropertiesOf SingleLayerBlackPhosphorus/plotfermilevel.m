clear all;
close all;
T=[];N=[];Ef=[];
T0=1;N0=1e11;
X1=-5;X2=5;X0=0;
for T0=5:5:300
    X0=fermilevel(T0,N0,X1,X2);
    T=[T T0];
    Ef=[Ef X0];
    fprintf(1,'T= %d N = %d Ef = %d \n',T0,N0,X0);
end
plot(T,Ef,'ob','linewidth',3);