clear all
close all

%% semidiscrtizzo il problema nello spazio con differenze finite di ordine 2

m=1800;
a=0;
b=pi/2;
x=linspace(a,b,m)';
h=(b-a)/(m-1);
A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));

%% condizioni al bordo
d=1/2;
A(1,1:2)=[0,0];
A(m,m-1:m)=[2,-2]/h^2;

%termine sorgente
b=@(t) (sin(x)/2)*( cos(t/2) - sin(t/2) ); %% non serve azzerare la prima componente, poiché la prima componente del vettore
% x, che è 0, annulla sin(x).

%% Risoluzione problema differenziale mediante metodo dei trapezi
%% La soluzione analitica è u(t,x)=sin(x)*cos(t/2)

y0=sin(x); %dato iniziale
count=0;
tsrange=10:10:100;
tstar=1;

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    t=0;
    y0=sin(x);
    y=y0;
    
    for n=1:ts
        y=(speye(m) - k*0.5*d*A)\(y + k*0.5*(d*A*y + b(t) + b(t+k) ) ); % non è necessario newton dal momento che il problema è lineare
        t=t+k;
    end
    
    err(count)=norm(y - (sin(x)*cos(t/2)),inf);
    
end

figure
plot(x,y,'*',x,sin(x)*cos(0.5),'r')
title('Soluzioni per t*=1')
legend('soluzione numerica','soluzione analitica')

figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'r')
title('Ordine di convergenza temporale')
legend('errore per t*=1','ordine 2')
xlabel('timesteps')
ylabel('err')
