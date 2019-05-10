clear all
close all

mrif=2^11+1;
m=mrif;
h=2/(m-1);
x=linspace(-1,1,m)';

D2 = toeplitz(sparse([1,1],[1,2],[-2,1]/(h^2),1,m));
D1 = toeplitz(sparse(1,2,-1/(2*h),1,m),sparse(1,2,1/(2*h),1,m));

A=diag(exp(-x))-D2;
b=exp(-x.^2);
A(1,1:2)=[1,0];
b(1)=0;
A(m,m-1:m)=[-2/h^2,2/h^2+exp(-1)];

u=A\b;

Urif=u;
count=0;
mrange=2.^(4:9)+1;

for m=mrange
  count++;
  h=2/(m-1);
  x=linspace(-1,1,m)';
  
  D2 = toeplitz(sparse([1,1],[1,2],[-2,1]/(h^2),1,m));
  D1 = toeplitz(sparse(1,2,-1/(2*h),1,m),sparse(1,2,1/(2*h),1,m));
  
  A=diag(exp(-x))-D2;
  b=exp(-x.^2);
  A(1,1:2)=[1,0];
  b(1)=0;
  A(m,m-1:m)=[-2/h^2,2/h^2+exp(-1)];
  
  u=A\b;
  
  err(count)=norm(u-Urif(1:(mrif-1)/(m-1):mrif),inf);
end

loglog(mrange,err,'*',mrange,err(end)*(mrange/mrange(end)).^(-2),'g')
