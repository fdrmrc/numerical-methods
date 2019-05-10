clear all
close all

f=@(y) [-sin(y(1))*y(2);sin(y(1))*y(1)];
df=@(y)[-cos(y(1))*y(2),-sin(y(1));cos(y(1))*y(1)+sin(y(1)),0];
ts=1000;
k=1/ts;
y=NaN(2,ts+1);
y0=[1;1];
y(:,1)=y0;
tol=k^2/100;
SC(1)=y0'*y0;
for n=1:ts
  F=@(x) x- y(:,n) -k*f(0.5*(y(:,n)+x));
  J=@(x) eye(2) - k*0.5*df(0.5*(y(:,n)+x));
  x0=y(:,n);
  res=-J(x0)\F(x0);
  while (norm(res,inf)>tol)
    x0+=res;
    res=-J(x0)\F(x0);
  end
  x0+=res;
  y(:,n+1)=x0;
  SC(n+1)=dot(y(:,n),y(:,n)); %scalar product. See documentation

end
t=linspace(0,1,ts+1);
figure
plot(t,SC,'*')
title('Andamento quantita nel tempo')
%E'costantemente 2, per una proprietà del metodo scelto. Se per il problema differenziale vale che y(t)^T y(t) è costante, allora vale che la quantità y_{n}^T y_{n} è costante 
