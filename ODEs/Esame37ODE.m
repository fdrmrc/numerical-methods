clear all
close all

f=@(y) sqrt((2-y)/y);
df=@(y) -1/(sqrt(2/y - 1)*(y^2));
%Fondamentale usare un metodo implicito, in modo da evitare la valutazione della f in y0=0
%Fondamentale anche la scelta del guess iniziale in Newton
tstar=pi/2;
tsrif=4000;
ts=tsrif;
k=tstar/ts;
y=NaN(1,ts+1);
y0=0;
y(1,1)=y0;
tol=k^2;
for n=1:ts
  F=@(x) x- y(1,n) -k*f(x);
  J=@(x) 1-k*df(x);
  x0=y(1,n)+k;
  res=-J(x0)\F(x0);
  while (norm(res,inf)>tol)
    x0+=res;
    res=-J(x0)\F(x0);
  end
  x0+=res;
  y(1,n+1)=x0;
end

yRif=y;
count=0;
tsrange=[10:10:100];

for ts=tsrange
  count++;
  k=tstar/ts;
  y=NaN(1,ts+1);
  y0=0;
  y(1,1)=y0;
  tol=k^2;
  for n=1:ts
    F=@(x) x- y(1,n) -k*f(x);
    J=@(x) 1-k*df(x);
    x0=y(1,n)+k;
    res=-J(x0)\F(x0);
    while (norm(res,inf)>tol)
      x0+=res;
      res=-J(x0)\F(x0);
    end
    x0+=res;
    y(1,n+1)=x0;
  end

  err(count)=norm(y(1,ts+1)-yRif(1,end),Inf);
endfor
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-1),'g')


