clear all; close all;

%% Risoluzione eseame 61.
%% E' un metodo Runge-Kutta semiimplicito, , con xi1=yn, a 3 stadi

%% definizione parametri
gamma=2-sqrt(2);
d=gamma/2;
omega=(1-d)/2;

%a
a(1,1)=0;
a(2,1)=d;
a(2,2)=d;
a(3,1)=omega;
a(3,2)=omega;
a(3,3)=d;
%b
b(1)=omega;
b(2)=omega;
b(3)=d;
%c
c(1)=0;
c(2)=gamma;
c(3)=1;

%% Dati del probl.
y0=0;
tsrif=1500;
ts=tsrif;
tstar=1;
k=tstar/ts;
y(1,1)=y0;
f=@(t,y) y+t; 


t=linspace(0,tstar,ts+1);
xi=zeros(1,3); % ho tre stadi
for n=1:ts
    
   xi(1,1)=y(1,n); %poiché a(1,1)=0=> non va risolto con Newton.
   
   xi(1,2)=(y(1,n) +d*k*(k*c(1)+k*c(2)+y(1,n)+2*t(n)))/(1-d*k); %jacobiano (cioè derivata) è una funzione lineare => NON serve newton
   
   xi(1,3)=( y(1,n) + omega*k*(2*t(n)+k*c(1)+xi(1,1) +k*c(2) + xi(1,2))+d*k*t(n)+d*k^2*c(3) )/(1-d*k);
    
   y(1,n+1)= y(1,n) + k*( b(1)*f(t(n)+k*c(1),xi(1,1)) + b(2)*f(t(n)+k*c(2),xi(1,2)) + b(3)*f(t(n)+k*c(3),xi(1,3)));
    
    
end
anal=@(x) -1 + exp(x) - x;
% 
% plot(t,y,'*',t,anal(t),'r')

%% Mostro convergenza
tsrange=[50:50:600];
counter=0;
for ts=tsrange
    counter=counter+1;
    k=tstar/ts;
    t=linspace(0,tstar,ts+1);
    y=NaN(1,ts+1);
    y(1,1)=y0;
    for n=1:ts
    
   xi(1,1)=y(1,n); %poiché a(1,1)=0=> non va risolto con Newton.
   
   xi(1,2)=(y(1,n) +d*k*(k*c(1)+k*c(2)+y(1,n)+2*t(n)))/(1-d*k); %jacobiano (cioè derivata) è una funzione lineare => NON serve newton
   
   xi(1,3)=( y(1,n) + omega*k*(2*t(n)+k*c(1)+xi(1,1) +k*c(2) + xi(1,2))+d*k*t(n)+d*k^2*c(3) )/(1-d*k);
    
   y(1,n+1)= y(1,n) + k*( b(1)*f(t(n)+k*c(1),xi(1,1)) + b(2)*f(t(n)+k*c(2),xi(1,2)) + b(3)*f(t(n)+k*c(3),xi(1,3)));
    
    
    end
    err(counter)=abs(y(1,ts+1)-anal(1));
end
figure
plot(t,y,'*',t,anal(t),'k')
legend('sol numerica','sol analitica')
figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'r')
