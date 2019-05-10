clear all
close all

mrif=2^11+1;
m=mrif;
h=2/(m-1);
x=linspace(-1,1,m)';

D2 = toeplitz(sparse([1,1],[1,2],[-2,1]/(h^2),1,m));
D1 = toeplitz(sparse(1,2,-1/(2*h),1,m),sparse(1,2,1/(2*h),1,m));

F=@(u) [-2*(2*u(2)-2*u(1))/h^2-u(1)*(1-u(1));(-2*(D2*u)+diag(x)*(D1*u)-u.*(1-u))(2:m-1);u(m)-2];
J=@(u) [[4/h^2-1+2*u(1),-4/h^2,zeros(1,m-2)];(-2*D2+diag(x)*D1-speye(m)+2*diag(u))(2:m-1,1:end);[zeros(1,m-1),1]];
u0=ones(m,1);
res=-J(u0)\F(u0);
tol=h^2;
while(norm(res,inf)>tol)
  u0+=res;
  res=-J(u0)\F(u0);
end
u0+=res;

Urif=u0;

mrange=2.^(4:9)+1;
count=0;
for m=mrange
  count++;
  h=2/(m-1);
  x=linspace(-1,1,m)';

  D2 = toeplitz(sparse([1,1],[1,2],[-2,1]/(h^2),1,m));
  D1 = toeplitz(sparse(1,2,-1/(2*h),1,m),sparse(1,2,1/(2*h),1,m));
  
  F=@(u) [-2*(2*u(2)-2*u(1))/h^2-u(1)*(1-u(1));(-2*(D2*u)+diag(x)*(D1*u)-u.*(1-u))(2:m-1);u(m)-2];
  J=@(u) [[4/h^2-1+2*u(1),-4/h^2,zeros(1,m-2)];(-2*D2+diag(x)*D1-speye(m)+2*diag(u))(2:m-1,1:end);[zeros(1,m-1),1]];
  u0=ones(m,1);
  res=-J(u0)\F(u0);
  tol=h^2;
  while(norm(res,inf)>tol)
    u0+=res;
    res=-J(u0)\F(u0);
  end
  u0+=res;

  err(count)=norm(u0-Urif(1:(mrif-1)/(m-1):mrif),inf);
end

loglog(mrange,err,'*',mrange,err(end)*(mrange/mrange(end)).^(-2),'g')

%Se u(1)=1, allora la soluzione sarebbe u(x)=1 \forall x \in (-1,1), da cui l'errore di discretizzazione è esattamente nullo. Quindi non è possibile mostrare l'ordine di convergenza
