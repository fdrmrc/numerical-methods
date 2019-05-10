clear all
close all

find=false;
mrange=2.^(4:9)+1;

for m=mrange  
  h=1/(m-1);
  x=linspace(1,2,m)';
  D1= toeplitz(sparse(1,2,-1/(2*h),1,m), sparse(1,2,1/(2*h),1,m));
  D2 = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m)); 
  A=D2-diag(2./x)*D1+diag(2./(x.^2)).*speye(m);
  A(1,1:2)=[1,0];
  A(m,m-1:m)=[0,1];
  b=sin(log(x));
  b(1)=1;
  b(m)=2;

  u=A\b;
  
  uPrimoDue=(3*u(m)-4*u(m-1)+u(m-2))/(2*h); %nota simmetria rispetto al calcolo di u'(x_1)

  tol=1e-4;
  if find==false && uPrimoDue<tol
    trovato =true
  endif

end

s=@(x) 1/2 *x.*(-x -  x.*sin(log(x)) + 2*x.* sin(log(2)) + x.*(-cos(log(x)))+... %soluzione analitica
		+ 2*x.*cos(log(2)) + 4 - 2*sin(log(2)) - 2*cos(log(2)));

plot(x,u,'*', x,s(x),'r',x,uPrimoDue*x -2*uPrimoDue +2,'g') %la seconda parte rappresenta una retta che ha la pendenza richiesta al bordo
xlabel('x')
ylabel('u(x)')

    
  
