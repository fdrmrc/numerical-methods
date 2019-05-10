clear all
close all

mrif=201;
m=mrif;
h=2/(m-1);
x=linspace(-1,1,m)';

D2 = toeplitz(sparse([1,1],[1,2],[-2,1]/(h^2),1,m));
D1 = toeplitz(sparse(1,2,-1/(2*h),1,m),sparse(1,2,1/(2*h),1,m));

F=@(u) [-2*(2*u(2)-2*u(1))/h^2 - sin(u(1));(-2*(D2*u)+diag(1-x)*(D1*u)-sin(u))(2:m-1);-2*(2*u(m-1)-2*u(m))/h^2 - sin(u(m))];
J=@(u) [[-cos(u(1))+4/h^2,-4/h^2,zeros(1,m-2)];(-2*D2+diag(1-x)*D1 - diag(cos(u)))(2:m-1,1:m);[zeros(1,m-2),-4/h^2,4/h^2-cos(u(m))]];

u0=pi*ones(m,1);

res=-J(u0)\F(u0);
tol=h^2;
while(norm(res,inf)>tol)
  u0+=res;
  res=-J(u0)\F(u0);
end
u0+=res;

plot(x,u0,'b-o')


%NON E' POSSIBILE MOSTRARE IL CORRETTO ORDINE SPAZIALE PERCHE' LA SOLUZIONE NON E' UNICA. Infatti u(x)=k*pi, k \in Z è una famiglia di soluzioni per l'equazione differenziale, e, essendo u(x) costante, le differenze finite sono esatte
