clear all
close all


m=101;
x=linspace(0,1,m)';
h=1/(m-1);

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
B = toeplitz(sparse(1,2,-1/(2*h),1,m),...
    sparse(1,2,1/(2*h),1,m));

A(1,1:2)=[-2,2]/h^2;
B(1,2)=0;
A(m,m-1:m)=[0,0];
B(m,m-1)=0;
d=2;

%termine reazione 
g=@(u) [u(1:m-1).*(1-u(1:m-1));0]; %azzero l'ultima componente per imporre la condizione di compatibilità
dg=@(u) [-2*u(1:m-1) + 1;0];

%% risoluzione sistema differenziale
y0=(x.^2).*(1-x);
tstar=0.5;
tsrif=1000;
ts=tsrif;
k=tstar/ts;
y=y0;
for n=1:ts
    Jn=d*A + spdiags(dg(y),0,m,m);
    P=phi1m(k*Jn);
    y=y+k*P*(d*A*y + g(y));
%     plot(x,y)
%     xlabel('x')
%     ylabel('u(t*,x)')
%     title('soluzione numerica')
%     pause(0.01)    
    
end

%% Mostro ordine di convergenza temporale

yRef=y;
count=0;
tsrange=10:10:100;

for ts=tsrange
    k=tstar/ts;
    count=count+1;
    y=y0;
for n=1:ts
    Jn=d*A + spdiags(dg(y),0,m,m);
    P=phi1m(k*Jn);
    y=y+k*P*(d*A*y + g(y));    
end

err(count)=norm(y-yRef,inf);

end

figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'g')
title('Error plot')
legend('error','ord2')