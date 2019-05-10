clear all
close all

m=101;
x=linspace(0,1,m)';
h=1/(m-1);
d=1/20;

A = toeplitz(sparse([1,1],[1,2],[-2,1]/(h^2),1,m));

A(1,1:2)=[0,0];
A(m,m-1:m)=[2,-2]/h^2;
% poiché x(1)=0, la prima riga è azzerata automaticamente. Pertanto
%non serve inserire una matrice identità e azzerarne la prima componente
%Dalle condizioni di compatibilità ricavo che u(t,0)=0.

%% Risoluzione problema differenziale
y0=x.*(1-x).^2;
y=y0;
tstar=1;
tsrif=500;
ts=tsrif;
k=tstar/ts;

g=@(u) x./(1+u.^2);
t=0;
for n=1:ts
    t=t+k;
    Jn=d*A + spdiags((-2*x.*y)./(1+y.^2).^2,0,m,m);
    P=phi1m(k*Jn);
    y=y+k*P*(d*A*y + g(y));
    plot(x,y)
    title(sprintf('Soluzione a t = %0.2f',t));
    axis([0,1,0,1])
    pause(0.001)
end
yRef=y;


tsrange=10:10:100;
count=0;

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    y0=x.*(1-x).^2;
    y=y0;
    
    g=@(u) x./(1+u.^2);
    dg=@(u) -(2*x.*u)./((1+u.^2).^2);
    
    for n=1:ts
        Jn=d*A + spdiags((-2*x.*y)./(1+y.^2).^2,0,m,m);
        P=phi1m(k*Jn);
        y=y+k*P*(d*A*y + g(y));
    end
    
    err(count)=norm(y-yRef,inf);
    
end

figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'g')
legend('err','ord2')