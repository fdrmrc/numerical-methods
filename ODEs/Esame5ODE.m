clear all
close all

f=@(t,y) y*(1-y);
df=@(t,y) -2*y+1;
y0=0.5;

tstar=1;
tsrif=1000;
ts=tsrif;
k=tstar/ts;
t=linspace(0,tstar,ts+1);
y=NaN(1,ts+1);
y(1,:)=y0;

for n=1:ts
  F=@(x) x-y(:,n)-k*f(t(n)+0.5*k,0.5*(y(:,n)+x));
  J=@(x) 1-k*0.5*df(t(n)+0.5*k,0.5*(y(:,n)+x));

  x0=y(:,n); %guess iniziale

  res=-J(x0)\F(x0);
  tol=0.01*k^2;
  while(norm(res,inf)>tol)
    x0+=res;
    res=-J(x0)\F(x0);
  endwhile
  x0+=res;

  y(:,n+1)=x0;

endfor
yRif=y; %soluzione di riferimento
tsrange=10:10:200;
count=0;

for ts=tsrange
  count++;  
  k=tstar/ts;
  t=linspace(0,tstar,ts+1);
  y=NaN(1,ts+1);
  y(1,:)=y0;
  
  for n=1:ts
    F=@(x) x-y(:,n)-k*f(t(n)+0.5*k,0.5*(y(:,n)+x));
    J=@(x) 1-k*0.5*df(t(n)+0.5*k,0.5*(y(:,n)+x));
    
    x0=y(:,n); %guess iniziale

    res=-J(x0)\F(x0);
    tol=0.01*k^2;
    while(norm(res,inf)>tol)
      x0+=res;
      res=-J(x0)\F(x0);
    endwhile
    x0+=res;
    
    y(:,n+1)=x0;
  
  endfor
  err(count)=norm(y(:,n+1)-yRif(:,end));
endfor

loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'g')



