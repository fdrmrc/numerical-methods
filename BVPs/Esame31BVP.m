clear all
close all

m=200;
h=1/(m-1);
x=linspace(0,1,m)';

D1 = toeplitz(sparse([1,1],[1,2],[-2,1]/(h^2),1,m)); 
D2= toeplitz(sparse(1,2,-1/(2*h),1,m), sparse(1,2,1/(2*h),1,m));

F=@(u) [u(1)^2-2+(-2*u(1)+2*u(2)-2*h)/(h^2);(D2*u-D1*u+u.^2-exp(2*x))(2:m-1);u(m)-exp(1)];
J=@(u) [[2*u(1)-2/h^2,2/h^2,zeros(1,m-2)];(D2-D1+2*speye(m))(2:m-1,1:m);[zeros(1,m-1),1]];
u0=ones(m,1);
tol=h^2;
res=-J(u0)\F(u0);
while (norm(res,inf)>tol)
  u0+=res;
  res=-J(u0)\F(u0);
end
u0+=res;

				%Verifica dell'approssimazione
mrange=2.^(4:10)+1;
count=0;
for m=mrange
  count++;
  h=1/(m-1);
  x=linspace(0,1,m)';

  D1 = toeplitz(sparse([1,1],[1,2],[-2,1]/(h^2),1,m)); 
  D2= toeplitz(sparse(1,2,-1/(2*h),1,m), sparse(1,2,1/(2*h),1,m));
  
  F=@(u) [u(1)^2-2+(-2*u(1)+2*u(2)-2*h)/(h^2);(D2*u-D1*u+u.^2-exp(2*x))(2:m-1);u(m)-exp(1)];
  J=@(u) [[2*u(1)-2/h^2,2/h^2,zeros(1,m-2)];(D2-D1+2*speye(m))(2:m-1,1:m);[zeros(1,m-1),1]];
  u0=ones(m,1);
  tol=h^2;
  res=-J(u0)\F(u0);
  while (norm(res,inf)>tol)
    u0+=res;
    res=-J(u0)\F(u0);
  end
  u0+=res;

  %CN=(u0(2)-u0(1))/h; %approx del prim'ordine
  CN=(-3*u0(1)+4*u0(2)-u0(3))/(2*h);
  err(count)=abs(CN-1);

endfor
loglog(mrange,err,'*',mrange,err(end)*(mrange/mrange(end)).^(-2),'g')
legend('approx. u''(0)=0','ord 2')
