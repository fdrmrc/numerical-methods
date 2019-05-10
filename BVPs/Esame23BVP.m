clear all
close all

mrif=2^12+1;
m=mrif;
h=1/(m-1);
x=linspace(1,2,m)';
D2 = toeplitz(sparse([1,2],[1,1],[-2,1]/(h^2),m,1));

F=@(u) [u(1)-2;(D2*u-2*u.^3+6*u+2*x.^3)(2:m-1);u(m)-2.5];
J=@(u) [[1,zeros(1,m-1)];(D2-6*diag(u.^2)+6*speye(m))(2:m-1,1:m);[zeros(1,m-1),1]];

u0=linspace(2,2.5,m)';
res=-J(u0)\F(u0);
tol=h^2;
while(norm(res,inf)>tol)
  u0+=res;
  res=-J(u0)\F(u0);
end
u0+=res;

Urif=u0;
count=0;
mrange=2.^(4:8)+1;

for m=mrange
  count++;
  h=1/(m-1);
  x=linspace(1,2,m)';
  D2 = toeplitz(sparse([1,2],[1,1],[-2,1]/(h^2),m,1));
  
  F=@(u) [u(1)-2;(D2*u-2*u.^3+6*u+2*x.^3)(2:m-1);u(m)-2.5];
  J=@(u) [[1,zeros(1,m-1)];(D2-6*diag(u.^2)+6*speye(m))(2:m-1,1:m);[zeros(1,m-1),1]];

  u0=linspace(2,2.5,m)';
  res=-J(u0)\F(u0);
  tol=h^2;
  while(norm(res,inf)>tol)
    u0+=res;
    res=-J(u0)\F(u0);
  end
  u0+=res;

  err(count)=norm(u0-Urif(1:(mrif-1)/(m-1):mrif),inf);
  
endfor

loglog(mrange,err,'o',mrange,mrange.^-2*(err(end)/mrange(end)^-2))
legend('errore','ord2')
