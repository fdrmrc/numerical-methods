clear all; close all;

%% Risoluzione problema diff non omogeneo con esponenziale punto medio
m=400;
h=pi/(m-1);
x=linspace(pi,2*pi,m)';

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));

%% condizioni al bordo
A(m,m-1:m)=[2,-2]/h^2;
A(1,1:2)=[0,0];

b=@(t) 5/4*exp(t)*cos(x/2); %non c'è bisogno di mettere in prima componente 0

%% convergenza
y0=cos(x/2);

counter=0;
tsrange=10:10:100;

for ts=tsrange
    counter=counter+1;
    k=1/ts;
    P=phi1m(k*A);
    y=y0;
    
    t=0;
    for n=1:ts
        y=y+k*P*(A*y+b(t+k/2));
        t=t+k;
    end
    
    err(counter)=norm(y-exp(t)*cos(x/2));
    
end

loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'r')

legend('errore','ordine2')
title('Error plot')
xlabel('time steps')
ylabel('error')
