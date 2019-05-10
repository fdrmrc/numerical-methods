clear all; close all;

%% Risoluzione sistema differenziale 79.

%% Runge-Kutta semi-implicito a tre stadi, ma con a(1,1)=0=> non è richiesto l'utilizzo di metodi numerici per la de
%% =terminazione dei coefficienti xi, che quindi possono essere ricavati risolvendo sistemi lineari 2x2.

%% definizione dei parametri:
theta=2-sqrt(2);
c(1)=0;
c(2)=theta;
c(3)=1;
    
b(1)=(3*theta-theta^2-1)/(2*theta);
b(2)=(1-theta)/(2*theta);
b(3)=theta/2;

%% dati del problema
ts=200;
tstar=2*pi;
k=tstar/ts;
A=[-500,0;0,-1]%matrice scalare che moltiplica il vttore incongnite.
r=@(t) [500*cos(t) - sin(t);sin(t)+cos(t)];%termine noto, dipendente dal tempo
y0=[1,0]';
f=@(t,y) A*y + r(t);

y=NaN(length(y0),ts+1);
y(:,1)=y0;

t=linspace(0,tstar,ts+1);
xi=zeros(length(y0),3); %matrice 2x3, ho tre stadi e il problema è bidimensionale

for n=1:ts
    xi(:,1)=y(:,n);
    xi(:,2)=(eye(2) - k*(theta/2)*A)\( y(:,n) + k*(theta/2)*(A*y(:,n) + r(t(n)+k*c(1))) + k*(theta/2)*r(t(n)+k*c(2)) );
    xi(:,3)=(eye(2) - b(3)*k*A)\( y(:,n) +...
        k*(b(1)*(A*y(:,n) + r(t(n)+k*c(1))) + b(2)*(A*xi(:,2)+r(t(n)+k*c(2)))) ...
        + b(3)*k*r(t(n)+k*c(3)));
    
    y(:,n+1)=y(:,n) + k*(b(1)*f(t(n)+k*c(1),xi(:,1)) + b(2)*f(t(n)+k*c(2),xi(:,2)) + b(3)*f(t(n)+k*c(3),xi(:,3)));
    
end

plot(t,y(1,:),'*',t,cos(t),'r',t,y(2,:),'sk',t,sin(t),'b')
legend('soluzione numerica 1','sol anal 1','sol numerica 2','sol anal 2')
figure
plot(y(1,:),y(2,:),'y-*') %nel piano della fasi ottengo una circonferenza