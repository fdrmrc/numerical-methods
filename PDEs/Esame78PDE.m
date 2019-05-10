clear all
close all


m=100;
h=1/(m-1);
x=linspace(0,1,m)';

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
B = toeplitz(sparse(1,2,-1/(2*h),1,m), sparse(1,2,1/(2*h),1,m));
d=1/10;
c=-10;

Pe=(abs(c)*h)/(2*d);
Pe<1

% condizioni ai bordi
A(1,1:2)=[-2,2]/h^2;
B(1,2)=0;
A(m,m-1:m)=[0,0];
B(m,m-1)=0;

%% Risoluzione mediante punto medio implicito
C=d*A-c*B;
y0=(x.^2).*(1-x);
tstar=0.05;
ts=1000;
k=tstar/ts;
y=y0;
t=0;
for n=1:ts
    t=t+k;
    y=(speye(m)-k*0.5*C)\((speye(m)+k*0.5*C)*y);
    plot(x,y,'r',x,x-x,'--')
    title(sprintf('Soluzione a t = %0.2f',t));
    axis([0,1,-0.3,1]);
    pause(0.01)
    for i=1:m
        if(y(i)<0)
        disp('pechlet negativo')
        end
    end
end
yrif=y;

%% verifica della convergenza
count=0;
tsrange=10:10:100;

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    y=y0;
    
    for n=1:ts
        y=(speye(m)-k*0.5*C)\((speye(m)+k*0.5*C)*y);
    end
    err(count)=norm(y-yrif,inf);
end

figure
loglog(tsrange,err,'*',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-1),'r',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'g')
legend('error','ord1','ord2')
