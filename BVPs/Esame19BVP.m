clear all
close all

m=200;
h=1/(m-1);
x=linspace(0,1,m)';

D2 = toeplitz(sparse([1,2],[1,1],[-2,1]/(h^2),m,1));

F=@(u) [2/(h^2)*(u(2)-u(1))-1-exp(u(1));(D2*u-exp(u)-2+exp(x.^2))(2:m-1);u(m)-1];
J=@(u) [[-2/(h^2)-exp(u(1)),2/(h^2),zeros(1,m-2)];(D2-diag(exp(u)))(2:m-1,1:end);[zeros(1,m-1),1]];

tol=0.01*h^2;
u0=ones(m,1);
res=-J(u0)\F(u0);
while(norm(res,Inf)>tol)
  u0+=res;
  res=-J(u0)\F(u0);
endwhile
u0+=res;

plot(x,u0,'r-o')

%La soluzione analitica è u(x)=x^2. Pertanto il metodo delle differenze finite è esatto (u''''(x)=0 \forall x \in (0,1) ), da cui l'impossibilità di mostare un grafico loglog dell'errore
