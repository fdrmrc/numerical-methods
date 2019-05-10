clear all
close all

count=0;
mrange=2.^(4:10)+1;

for m=mrange
  count++;
  h=(2*pi)/(m-1);
  x=linspace(0,2*pi,m)';
  
  D1= toeplitz(sparse(1,2,-1/(2*h),1,m), sparse(1,2,1/(2*h),1,m));
  D2 = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
  
  A=D2+3*D1-speye(m);
  b=3*cos(x)-2*sin(x);
  
  A(1,1:2)=[-1-2/h^2,2/h^2];
  b(1)=+2/h;
  A(m,m-1:m)=[0,1];
  b(m)=0;
  u0=A\(b);

  err(count)=norm(u0-sin(x),Inf);
endfor

loglog(mrange,err,'*',mrange,err(end)*(mrange/mrange(end)).^(-2),'g')
