clear all
close all

m=400;
a=0;
b=pi/2;
h=(b-a)/(m-1);
x=linspace(a,b,m)';
d=0.5;


A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
B = toeplitz(sparse(1,2,-1/(2*h),1,m),sparse(1,2,1/(2*h),1,m));

A(1,1:2)=[0,0];
B(1,2)=0;
A(m,m-1:m)=[2,-2]/h^2;
B(m,m-1)=0;

b=@(t) [0;exp(-t/2)*cos(x(2:m))];
C=d*A-B;
count=0;
tstar=1;
tsrange=10:10:100;
y0=sin(x);

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    P=phi1m(k*C);
    y=y0;
    t=0;
    for n=1:ts
        y=y + k*P*(C*y+b(t+k/2));
        t=t+k;
    end
    
    err(count)=norm(y-exp(-t/2)*sin(x),inf);
    
end
figure
plot(x,y)
figure
loglog(tsrange,err,'g*',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'r')


