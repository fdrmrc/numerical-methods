clear all
close all

m=300;
h=1/(m-1);
x=linspace(0,1,m)';

d=0.0025;
c=-1;
p=10;
A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
B = toeplitz(sparse(1,2,-1/(2*h),1,m), sparse(1,2,1/(2*h),1,m));

%% condizioni al bordo
A(1,1:2)=[-2,2]/h^2;
B(1,2)=0;

A(m,m-1:m)=[0,0];
B(m,m-1)=0;

%! va azzerato il termine di reazione

g=@(u) [p*(u(1:m-1).^2).*(1-u(1:m-1));0];
dg=@(u) [p*(2-3*u(1:m-1)).*u(1:m-1);0];

%% Risoluzione sistema differenziale

y0=x.^2.*(1-x);
y=y0;
tstar=0.5;
ts=800;
k=tstar/ts;
C=d*A-c*B;


for n=1:ts
   Jn=C + diag(dg(y));
   P=phi1m(k*Jn);
   y=y+k*P*(C*y + g(y)); 
   %plot(x,y,'r',x,y0,'o')
   %axis([0,0.1,0,1])
   %pause(0.01)
end

count=0;
tsrange=10:10:100;
yrif=y;

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    y=y0;
    for n=1:ts
        Jn=C + diag(dg(y));
        P=phi1m(k*Jn);
        y=y+k*P*(C*y + g(y));
    end
    err(count)=norm(y-yrif,inf);
    
end

figure
loglog(tsrange,err,'*',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-2))
title('error plot')
xlabel('timesteps')
ylabel('error')
legend('error','ord2')
