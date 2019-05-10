clear all
close all

mrif=1001;
m=mrif;
h=1/(m-1);
x=linspace(0,1,m)';

d=0.005;
c=3;
A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
B = toeplitz(sparse(1,2,-1/(2*h),1,m), sparse(1,2,1/(2*h),1,m));

%condizioni al bordo
A(1,1:2)=[0,0];
B(1,1:2)=[0,0];

A(m,m-1:m)=[2,-2]/h^2;
B(m,m-1:m)=[0,0];
Pe=(c*h)/(2*d);
Pe<1
 

C=d*A - c*B;


%risoluzione problema differenziale y'=C*y
y0=x.*(1-x).^2;
tsrif=200;
ts=tsrif;
tstar=1;

%trapezi

k=tstar/ts
RestrTrap=(h^2)/d
k<RestrTrap

yTrap=y0;
for n=1:ts
	yTrap=(speye(m) - k*0.5*C)\( (speye(m) + k*0.5*C )*yTrap);
    plot(x,yTrap)
    pause(0.005)
    axis([0,1.5,-1,1])
end

%a causa del passo e della restrizione sul passo di trapezi,
%questa soluzione non è fisicamnente accettabile (nonostante Pe<1)
    

%eulero implicito
yEI=y0;

for n=1:ts
	yEI=(speye(m)-k*C)\(yEI);
    plot(x,yEI)
   axis([0,1,0,0.16])
   pause(0.1)
end

%la soluzione con eulero implicito è fisicamente più accettabile, dal
%momento che non presenta oscillazioni e si mantiene positiva

yrif=yEI;
