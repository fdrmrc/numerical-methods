clear all
close all

				%metodo utilizzato: Eulero implicito
  %risolvo prima il secondo sistema per determinare numericamente y(t)

%parametri
b=0.01;
m=0.5;
g=9.81;
alpha=pi/3; %radianti 
v0=200;
B=b/m;

A=[0,1;0,-B];
b=[0;-g];
tstar=40;
ts=2000;
k=tstar/ts;
y0=[0;v0*sin(alpha)];
y=NaN(2,ts+1);
y(:,1)=y0;

for n=1:ts
  y(:,n+1)=(eye(2)-k*A)\(y(:,n)+k*b);
end
t=linspace(0,tstar,ts+1);
plot(t,y(1,:))

%ricerca zero di y(t)
prod=eps;
j=1;
while prod>=0
  prod=y(1,j)*y(1,j+1);
  j++;
endwhile

		      %scelgo come zero la media tra questi due valori
jstar=((j-1)+j)/2;
T=jstar*k; %Tstar=t0+j*k: approssimazione del tempo T tale che y(T)=0

%% Risoluzione secondo problema
x0=[0;v0*cos(alpha)];
AA=[0,1;0,-B];
x=NaN(2,ts+1);
x(:,1)=x0;

for n=1:ts
  x(:,n+1)=(eye(2)-k*AA)\(x(:,n));
end
x(1,floor(jstar)) %punto di atterraggio


  
