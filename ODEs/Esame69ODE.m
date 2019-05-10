clear all; close all;
%% Risoluzione esame 69 con Eulero-Rosenbrock esponenziale
c=-20;
tstar=1;
tsrif=5000;
ts=tsrif;
t=linspace(0,tstar,ts+1);
y0=[0,0]';
k=tstar/ts;

%% definizione funzione
f=@(y) [c*(y(1)-cos(y(2)));1];
anal=@(x) 20/401*(-20*exp(-20*x) + 20*cos(x) + sin(x));

y=NaN(2,ts+1);
y(:,1)=y0;

for n=1:ts
    Jn=[c,c*sin(y(2,n));0,0];
    y(:,n+1)=y(:,n) + k*phi1m(k*Jn)*f(y(:,n));
end
figure
plot(t,y(1,:),'b-s',t,anal(t),'y')
legend('soluzione numerica','sol analitica')
yrif=y;

%% ordine di convergenza rispetto a soluzione di riferimento

counter=0;
tsrange=[50:50:700];

for ts=tsrange
counter=counter+1;    
y0=[0,0]';
y=NaN(2,ts+1);
y(:,1)=y0;
k=tstar/ts;
t=linspace(0,tstar,ts+1);

for n=1:ts
    Jn=[c,c*sin(y(2,n));0,1];
    y(:,n+1)=y(:,n) + k*phi1m(k*Jn)*f(y(:,n));
end
err(counter)=norm(y(1,ts+1) -anal(tstar),inf); % yrif(1,tsrif+1)
errR(counter)=norm(y(1,ts+1) - yrif(1,tsrif+1),inf);
end


figure
loglog(tsrange,err,'g*',tsrange,err(end)*(tsrange/tsrange(end)).^(-1),'r',...
    tsrange,errR,'b*')

legend('errore rispetto sol anal','ordine 1','errore risp sol rif')
