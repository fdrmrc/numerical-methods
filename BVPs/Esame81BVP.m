clear all
close all

mrif=2^11+1;
m=mrif;
h=2/(m-1);
x=linspace(0,2,m)';

D2 = toeplitz(sparse([1,1],[1,2],[-2,1]/(h^2),1,m));
D1 = toeplitz(sparse(1,2,-1/(2*h),1,m),sparse(1,2,1/(2*h),1,m));

F=@(u) [u(1)-1;(-diag(x+1)*(D2*u)+diag(x.^2)*(D1*u)-1./(1+u.^2))(2:m-1);...
	-3*(-2*u(m)+2*u(m-1))/h^2-1./(1+u(m)^2)];
J=@(u) [[1,zeros(1,m-1)];...
	(-diag(x+1)*D2+diag(x.^2)*D1+diag(2*u./((1+u.^2).^2)))(2:m-1,1:m);...
	[zeros(1,m-2),-6/h^2,6/h^2+(2*u(m))/((1+u(m))^2)]];
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
  
  x=linspace(0,2,m)';
  
  D2 = toeplitz(sparse([1,1],[1,2],[-2,1]/(h^2),1,m));
  D1 = toeplitz(sparse(1,2,-1/(2*h),1,m),sparse(1,2,1/(2*h),1,m));
  
  F=@(u) [u(1)-1;(-diag(x+1)*(D2*u)+diag(x.^2)*(D1*u)-1./(1+u.^2))(2:m-1);...
	  -3*(-2*u(m)+2*u(m-1))/h^2-1./(1+u(m)^2)];
  J=@(u) [[1,zeros(1,m-1)];...
	  (-diag(x+1)*D2+diag(x.^2)*D1+diag(2*u./((1+u.^2).^2)))(2:m-1,1:m);...
	  [zeros(1,m-2),-6/h^2,6/h^2+(2*u(m))/((1+u(m))^2)]];
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
