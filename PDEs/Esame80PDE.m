clear all
close all


m=100;
h=1/(m-1);
x=linspace(0,1,m)';

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
B = toeplitz(sparse(1,2,-1/(2*h),1,m), sparse(1,2,1/(2*h),1,m));
d=1/10;
c=10;

Pe=(abs(c)*h)/(2*d);
Pe<1;

A(1,1:2)=[0,0];
B(1,2)=0;
A(m,m-1:m)=[2,-2]/h^2;
B(m,m-1)=0;

y0=x.*(1-x).^2;
tstar=0.05;
ts=1000;
k=tstar/ts;
y=y0;
C=d*A-c*B;
t=0;
for n=1:ts
    y=(speye(m)-k*C)\(y);
    t=t+k;
    plot(x,y)
    title(sprintf('Soluzione a t = %0.2f',t));
  	pause(0.01)
end

yrif=y;

count=0;
tsrange=10:10:100;

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    y=y0;
    
for n=1:ts
    y=(speye(m)-k*C)\(y);
end

err(count)=norm(y-yrif,inf);
end

figure
loglog(tsrange,err,'*',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-1),'r',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'g')
title('Errore')

xlabel('passi')
ylabel('errore')
legend('errore','ordine 1','ordine 2')