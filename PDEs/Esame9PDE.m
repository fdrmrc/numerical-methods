clear all; close all 

%% Risoluzione esame PDE tramite metodo delle linee

m=150;
h=1/(m-1);
x=linspace(0,1,m)';
A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));

%% condizioni ai bordi
A(1,1:2)=[0,0];
A(m,m-1:m)=[0,0];
b=@(t) [0;(2*x(2:m-1) - x(2:m-1).^2).*cos(t) + 2*sin(t);cos(t)];


%% convergenza temporale mediante eulero implicito
y0=zeros(m,1);
counter=0;
tsrange=[10:10:100];

for ts=tsrange
    counter=counter+1;
    tstar=1;
    k=tstar/ts;
    y=y0;
    t=0;
    for n=1:ts
        y=(eye(length(y0)) - k*A)\(y+k*b(t+k));
        t=t+k;
    end
    
    err(counter)=norm(y-(2*x - x.^2)*sin(1),inf); 
    
end
figure
plot(x,y)
figure
loglog(tsrange,err,'k*',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-1),'r'...
    ,tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'b')

legend('errore','ordine 1','ordine 2')