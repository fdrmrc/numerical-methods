clear all; close all

m=500;
h=(pi/2)/(m-1);
x=linspace(0,pi/2,m)';

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));

%% condizioni al bordo
A(1,1:2)=[0,0]; %dirichlet
A(m,m-1:m)=[2,-2]/h^2; %neumann
b=@(t) [2*exp(t)*sin(x)];


%% risoluzione problema differenziale
counter=0;
tsrange=10:10:100;
y0=sin(x);

for ts = tsrange
    counter = counter+1;
    k = 1/ts;
    y = y0;
    t = 0;
    P = phi1m(k*A);
    for n = 1:ts
      y = y+k*P*(A*y+b(t)); % eulero esponenziale
      t = t+k;
    end
    err(counter) = norm(y-exp(t)*sin(x),inf);
end

loglog(tsrange,err,'*',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-1))
    title('grafico dell''errore')
    
