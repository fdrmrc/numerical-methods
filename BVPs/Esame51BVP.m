clear all
close all

mrange=2.^(4:10)+1;
count=0;
for m=mrange
  count++;
  h=1/(m-1);
  x=linspace(0,1,m)';
  
  D1= toeplitz(sparse(1,2,-1/(2*h),1,m), sparse(1,2,1/(2*h),1,m));
  D2 = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
  
  A=D2+D1+diag(exp(-x));
  b=exp(x);
				%condizioni ai bordi
  A(1,1:2)=[(2*h-2)/h^2,2/h^2];
				%b(1)=1 è gia verificata
  A(m,m-1:m)=[0,1];
  b(m)=1;
  u=A\b;
  
				%Verifico condizione di Robin
  CN(count)=u(1)+(-3*u(1)+4*u(2)-u(3))/(2*h);
endfor

loglog(mrange,CN,'*',mrange,CN(end)*(mrange/mrange(end)).^(-2),'g')
  
