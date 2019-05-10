clear all; close all;

tstar=1;
tsrif=2000;
ts=tsrif;

k=tstar/ts;

f=@(y) [-2*y(1)*y(2);y(1)^2-y(2)^2+y(3)-1;-4*(y(1)+y(2))*y(3)];
J=@(y) [-2*y(2),-2*y(1),0;2*y(1),-2*y(2),+1;-4*y(3),-4*y(3),-4*(y(1)+y(2))];

y0=[0.5,2,sqrt(10)]';

y=NaN(length(y0),ts+1);
y(:,1)=y0;

for n=1:ts
y(:,n+1)=y(:,n) + k*phi1m(k*J(y(:,n)))*f(y(:,n));
end
t=linspace(0,tstar,ts+1);

plot(t,y)
yref=y;

%% CONVERGENZA
counter=0;
tsrange=[200:50:600];

for ts=tsrange
counter=counter+1;
k=tstar/ts;
y(:,1)=y0;

for n=1:ts
y(:,n+1)=y(:,n) + k*phi1m(k*J(y(:,n)))*f(y(:,n));
end

err(counter)=norm(yref(:,tsrif+1) - y(:,ts+1));

end

figure
loglog(tsrange,err,'*',...
tsrange,err(end)*(tsrange/tsrange(end)).^(-1),...
tsrange,err(end)*(tsrange/tsrange(end)).^(-2))
legend('Errore','Ordine 1','Ordine 2')
