clear all
close all

%u=y(1),v=y(2),w=y(3)
f=@(t,y) [-2*y(1)*y(2);y(1)^2+y(3)^2-y(2)^2-1;-2*y(3)*(y(1)+y(2))];
df=@(t,y) [-2*y(2),-2*y(1),0;2*y(1),-2*y(2),2*y(3);-2*y(3),-2*y(3),-2*(y(1)+y(2))];

y0=[1;2;15];
tsrif=1000;
ts=tsrif;
k=1/ts;
y=NaN(3,ts+1);
t=linspace(0,1,ts+1);
y(:,1)=y0;

for n=1:ts
  F=@(x) x-y(:,n)-k*0.5*(f(t(n),y(:,n))+f(t(n+1),x));
  JF=@(x) eye(3)-k*0.5*df(t(n+1),x);

  x0=y(:,n)+k*f(t(n),y(:,n)); %guess iniziale con Eulero esplicito
  tol=0.01*k^2; %trapezi ha errore globale prop a k^2
  res=-JF(x0)\F(x0);
  iter=0;
  while(norm(res,Inf)>tol)
    x0+=res;
    res=-JF(x0)\F(x0);
    iter++;
  endwhile
  x0+=res;

  y(:,n+1)=x0;
endfor

yRif=y;

tsrange=[50:10:200];
count=0;

for ts=tsrange
  count++;
  y=NaN(3,ts+1);
  t=linspace(0,1,ts+1);
  k=1/ts;
  y(:,1)=y0;
  
  for n=1:ts
    F=@(x) x-y(:,n) - k*0.5*(f(t(n),y(:,n))+f(t(n+1),x));
    JF=@(x) eye(3)-k*0.5*df(t(n+1),x);
    
    x0=y(:,n)+k*f(t(n),y(:,n)); %guess iniziale con Eulero esplicito
    tol=k^2/100; %trapezi ha errore globale prop a k^2
    res=-JF(x0)\F(x0);
    while(norm(res,Inf)>tol)
      x0+=res;
      res=-JF(x0)\F(x0);
    endwhile
    x0+=res;
    
    y(:,n+1)=x0;
  endfor
  err(count)=norm(y(:,ts+1)-yRif(:,tsrif+1),Inf);
endfor

loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'g')
title('Error plot')
xlabel('time steps')
ylabel('inf norm of error')

  
