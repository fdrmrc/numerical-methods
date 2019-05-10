clear all
close all

%% Risoluzione tramite eulero rosenbrock-exp. 
%% Il problema deve essere reso autonomo e riportato al prim'ordine

f=@(y) [y(1)+y(2) - cos(y(4));y(3);y(3)+y(1).*(y(2).^2);1]; %% y4'=1;

y0=[1,0,1,0]'; %y4(t0)=t0=0

tstar=1;
tsrif=2000;
ts=tsrif;
k=tstar/ts;
y=NaN(length(y0),ts+1);
y(:,1)=y0;


for n=1:ts
    Jn=@(y) [1,1,0,sin(y(4));...
            0,0,1,0;...
            y(2).^2,2*y(1).*y(2),1,0;...
            0,0,0,0];
    
    P=phi1m(k*Jn(y(:,n)));
    y(:,n+1)=y(:,n) + k*P*(f(y(:,n)));
end
t=linspace(0,tstar,ts+1);

figure
plot(t,y(1,:),'g',t,y(2,:),'r',t,y(3,:),'k',t,y(4,:),'b')
title('Rappresentazione delle varie soluzioni')
xlabel('t')
ylabel('y(t)')
legend('y(1)','y(2)','y(3)','y(4)')



%% Mostro l'ordine di convergenza aspettato, che è 2, per questo metodo. Scelgo come soluzione di riferimento quella ottenta
%% con 2000 passi.

yRef=y;
tsrange=200:100:600;
count=0;

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    y=NaN(length(y0),ts+1);
    y(:,1)=y0;


for n=1:ts
    Jn=@(y) [1,1,0,sin(y(4));...
            0,0,1,0;...
            y(2).^2,2*y(1).*y(2),1,0;...
            0,0,0,0];
    
    P=phi1m(k*Jn(y(:,n)));
    y(:,n+1)=y(:,n) + k*P*(f(y(:,n)));
end

err(count)=norm(yRef(:,tsrif+1)-y(:,ts+1),inf);

end

figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'g')
title('Errore in norma infinito')
legend('errore','ordine 2')
xlabel('timesteps')
ylabel('y')
