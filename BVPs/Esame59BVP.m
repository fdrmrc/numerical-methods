clear all
close all

mrange=2.^(4:10)+1;
count=0;
for m=mrange
  count++;
  h=2/(m-1);
  x=linspace(-1,1,m)';
  
  D2 = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
  D1 = toeplitz(sparse(1,2,-1/(2*h),1,m), sparse(1,2,1/(2*h),1,m));

  F=@(u)[u(1)-1;(D2*u+u.^2-15/4*sqrt(abs(x))-abs(x).^5)(2:m-1);u(m)-1];
  J=@(u)[[1,zeros(1,m-1)];(D2+2*diag(u))(2:m-1,1:m);[zeros(1,m-1),1]];
  u0=ones(m,1);

  tol=h^2;
  res=-J(u0)\F(u0);
  while(norm(res,inf)>tol)
    u0+=res;
    res=-J(u0)\F(u0);
  end
  u0+=res;

  err(count)=norm(u0-abs(x).^(5/2),inf);
endfor
loglog(mrange,err,'*',mrange,err(end)*(mrange/mrange(end)).^(-2),'g',...
       mrange,err(end)*(mrange/mrange(end)).^(-1),'r',mrange,err(end)*(mrange/mrange(end)).^(-1.5),'k')
legend('err','ord2','ord1','ord1.5')

%E' corretto che l'ordine di Convergenze non sia 2: infatti la soluzione non è nemmeno C1, in particolare, non è derivabile in x=0
