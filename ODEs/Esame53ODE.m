clear all; close all

%% Esame 53

%soluzione analitica 
anal=@(x) 20/401*(-20*exp(-20*x) + 20*cos(x) + sin(x));

%% parametri del tableaux
a(1,1)=(3+sqrt(3))/6;
a(2,1)=-sqrt(3)/3;
a(2,2)=(3+sqrt(3))/6;

c(1)=(3+sqrt(3))/6;
c(2)=(3-sqrt(3))/6;

b(1)=0.5;
b(2)=0.5;
%% dati del problema

tstar=1;
ts=150;
t=linspace(0,tstar,ts+1);
k=tstar/ts;
y=NaN(1,ts+1);
y0=0;
y(1,1)=y0;

f=@(t,y) -20*(y-cos(t)); %il coefficiente davanti impone di stare attenti ad un possibile problema stiff.
xi=zeros(1,2);
counter=0;
tsrange=[100:50:500];
for ts=tsrange
    counter=counter+1;
    t=linspace(0,tstar,ts+1);
    k=tstar/ts;
    y(1,1)=y0;
    
for n=1:ts
    xi(1,1)=(y(1,n) + 20*k*a(1,1)*cos(t(n)+k*c(1)))/(1+20*a(1,1)*k);
    xi(1,2)=(y(1,n) + k*a(2,2)*20*cos(t(n)+k*c(2)) + k*(a(2,1)*(-20*xi(1,1)+20*cos(t(n)+k*c(1)))))/(1+20*k*a(2,2));
    y(1,n+1)=y(1,n) + k*(b(1)*f(t(n)+k*c(1),xi(1,1)) + b(2)*f(t(n)+k*c(2),xi(1,2)));
    
end

    err(counter)=abs(anal(tstar)-y(1,ts+1));
    
end

plot(t,y,'*',t,anal(t),'r')
figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-3),'r',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'b')
legend('errore','ordine 3','ordine 2')
