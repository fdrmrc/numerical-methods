clear all
close all

				%Eulero implicito
f=@(t,y) [y(2);cos(y(1))-y(2)/t];
df=@(t,y) [0,1;-sin(y(1)),-1/t];
tsrif=3000;
ts=tsrif;
t=linspace(0,1,ts+1);
k=1/ts;
y0=[1;0];
% y0=[pi/2;0]; %! per seconda richiesta: decommentare
y=NaN(2,ts+1);
y(:,1)=y0;

for n=1:ts
  F=@(x) x-y(:,n)-k*f(t(n+1),x);
  J=@(x) eye(2)-k*df(t(n+1),x);
  x0=y(:,n);
  tol=k^2;
  res=-J(x0)\F(x0);
  while (norm(res,inf)>tol)
    x0+=res;
    res=-J(x0)\F(x0);
  endwhile
  x0+=res;
  y(:,n+1)=x0;
endfor
yRif=y;

count=0;
tsrange=[10:10:200];
for ts=tsrange
  count++;
  k=1/ts;
  y=NaN(2,ts+1);
  y(:,1)=y0;
  t=linspace(0,1,ts+1);
  for n=1:ts
    F=@(x) x-y(:,n)-k*f(t(n+1),x);
    J=@(x) eye(2)-k*df(t(n+1),x);
    x0=y(:,n);
    tol=k^2;
    res=-J(x0)\F(x0);
    while (norm(res,inf)>tol)
      x0+=res;
      res=-J(x0)\F(x0);
    endwhile
    x0+=res;
    y(:,n+1)=x0;
  endfor

  err(count)=norm(y(1,ts+1)-yRif(1,tsrif+1),Inf);

endfor

  
figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-1),'g')
title('Error plot')
legend('err','ordine 1')

%Nel caso in cui y(0)=pi/2, basta decommentare dove richiesto. Evidentemente le soluzioni sono y(t)=pi/2,y'(t)=0, per cui il metodo numerico è esatto e il grafico dell'errore non ha alcun senso
