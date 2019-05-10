clear all
close all

f=@(y) [y(2);-9000*y(2)-5000*(sin(y(1))+y(1))]; %definizione funzione
tstar=10;
ts=100000; %problema presenta aree di stiffness, per cui un metodo esplicito richiede una restrizione sul passo temporale piuttosto severa.
k=tstar/ts;
y=NaN(2,ts+1);
y0=[pi;0];
y(:,1)=y0;

for n=1:ts
  y(:,n+1)=y(:,n)+k*f(y(:,n));
end

t=linspace(0,10,ts+1);
plot(t,y,'-*')
title('Soluzioni con eulero esplicito')

%% (NON RICHIESTA) !!
%Risoluzione mediante uno schema implicito.

df=@(y) [0,1;-5000*(1+cos(y(1))),-9000];
ts=100;
k=tstar/ts;
y=NaN(2,ts+1);
y(:,1)=y0;
tol=k^2;

for n=1:ts
  F=@(x) x- y(:,n)-k*f(x);
  J=@(x) eye(2) -k*df(x);
  x0=y(:,n);
  res=-J(x0)\F(x0);
  while(norm(res,Inf)>tol)
    x0+=res;
    res=-J(x0)\F(x0);
  end
  x0+=res;
  y(:,n+1)=x0;
endfor
t=linspace(0,tstar,ts+1);
figure
plot(t,y,'-*')

title('Soluzioni con eulero implicito')
%%nota: con un passo molto più grande ho ottenuto una soluzione accettabile. Ovviamente, al bordo x=0, essendo il passo maggiore, non viene riprodotta come con eulero esplicito.
