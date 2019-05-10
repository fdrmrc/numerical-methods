clear all
close all

m=401;
x=linspace(0,pi,m)';
h=pi/(m-1);


A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
d=4;

%condizioni al bordo
A(1,1:2)=[-2,2]/h^2;
A(m,m-1:m)=[0,0];

b=@(t) 3*exp(2*t)*cos(x/2); %non è necessario azzerare ultima componente

%% Risoluzione sistema differenziale
tsrange=10:10:100;
tstar=1;
count=0;
y0=cos(x/2);

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    t=0;
    P=phi1m(k*d*A);
    y=y0;
    for n=1:ts
        y=y+k*P*(d*(A*y) + b(t+k/2));
        t=t+k;
    end
    plot(x,y)
    %axis([0,1])
    pause(0.1);
    
    err(count)=norm(y - exp(2*t)*cos(x/2),inf);
    
end

figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'g')
