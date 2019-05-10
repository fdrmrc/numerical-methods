clear all
close all

m=200;
h=2/(m-1);
x=linspace(1,3,m)';

D2 = toeplitz(sparse([1,2],[1,1],[-2,1]/(h^2),m,1));
D1 = toeplitz(sparse(2,1,-1/(2*h),m,1), sparse(2,1,1/(2*h),m,1));

				%imposizione condizioni ai bordi
			       %F=@(u) D2*u-(1/8)*(32+2*x^3-u.*(D1*u))
F=@(u) [u(1)-17;(D2*u-(1/8)*(32*ones(m,1)+2*x.^3-(D1*u).*u))(2:m-1);u(m)-43/3];
J=@(u) [[1,zeros(1,m-1)];(D2+(1/8)*(diag(D1*u)+diag(u)*D1))(2:m-1,1:m);[zeros(1,m-1),1]];

u0=linspace(17,43/3,m)';
J0=J(u0); %per Newton inesatto
res=-J(u0)\F(u0);
Tol=1e-12;

iter=0;
while(norm(res,inf)>Tol)
  u0+=res;
%  res=-J(u0)\F(u0);
  res=-J0\F(u0); %per Newton inesatto
  iter++;
endwhile

u0+=res;
plot(x,u0,'b-o',x,x.^2 +16./x,'r')
title('Soluzioni')
	

	
	
	
