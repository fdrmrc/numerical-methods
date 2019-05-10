clear all; close all;

%% Risoluzione esame PDE numero 13, con metodo delle linee.

%% Costruisco soluzione di riferimento.
epsilon=1/100;
rho=10; %parametri del problema

h=1/100; %dato dal testo
m=101; %h=1/(m-1)
x=linspace(0,1,m);
x=x';
Pe=h/(2*epsilon)
A = toeplitz(sparse([1,2],[1,1],[-2,1]/h^2,m,1));
B = toeplitz(sparse(2,1,-1/(2*h),m,1),...
    sparse(2,1,1/(2*h),m,1));

%% condizioni al bordo (Neumann omogenee)
%! in prima riga sono [-2,2]: il segno meno sempre agli "estremi"
A(1,1:2)=[-2,2]/h^2; 
A(m,m-1:m)=[2,-2]/h^2;

B(1,1:2) = [0,0];
B(m,m-1:m) = [0,0];

%% termine reazione
b=@(u) rho*u.*(u-0.5).*(1-u);

%% ora devo risolvere il problema differenziale y'=epsilon*A*y+B*y+reazione, con y(t0)=y0;

y0 = 10 * x.^2.*(1-x).^2 + 0.5;
k=1/100; %lo dava il testo
tsrif=100; %k=tstar/ts = 1/ts = 100
ts=tsrif;
C=epsilon*A + B; %la matrice del problema y'=Ay + b
P = phi1m(k*C);

y=y0;
for n=1:ts
    y = y + k*P*(C*y + b(y)); % eulero esponenziale
        plot(x,y)
        axis([0,1,0,1.5])
        pause(0.1)
end
yrif=y;

tsrange=5:5:25;
count=0;

for ts=tsrange
    count=count+1;
    y0 = 10 * x.^2.*(1-x).^2 + 0.5;
    k=1/ts; %lo dava il testo
    P = phi1m(k*C);
    
    y=y0;
    for n=1:ts
        y = y + k*P*(C*y + b(y)); % eulero esponenziale
    end
    
    
    err(count)=norm(y-yrif,inf);
    
end
figure


loglog(tsrange,err,'*',...,
    tsrange,err(end)*(tsrange/tsrange(end)).^(-1),'g',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'r')
legend('errore','ord1','ord2')
