clear all; close all;

%% Risoluzione esame 30 con FD + Eulero implicito

d=0.01;
c=8; %coefficienti
rho=50;

m=101;
h=1/(m-1);
x=linspace(0,1,m)';

Pe=(c*h)/(2*d);
Pe<1 

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
B = toeplitz(sparse(1,2,-1/(2*h),1,m),...
    sparse(1,2,1/(2*h),1,m));


%% condizioni al bordo
A(1,1:2)=[0,0];
B(1,2)=0;
A(m,m-1:m)=[2,-2]/h^2;
B(m,m-1)=0;
C=d*A - c*B;

%termine reazione
g = @(u) rho*u.*(u-1/2).*(1-u);
dg = @(u) rho*((u-1/2).*(1-u)+u.*(1-u)-u.*(u-1/2));


%% Risoluzione sistema differenziale
y0=(x-1).^2;
ts=100;
tstar=0.1;
k=tstar/ts;

F = @(y,yn) y-yn-k*C*y-k*[0;g(y(2:m))]; % Neumann
J = @(y) eye(length(y))-k*C-k*spdiags([0;dg(y(2:m))],0,m,m); % Neumann


yn=y0;
tol=k^2/100;

for n=1:ts
    y=yn;
    delta=-J(y)\F(y,yn);
    while(norm(delta,inf)>tol)
        y=y+delta;
        delta=-J(y)\F(y,yn);
    end
    y=y+delta;
    yn=y;
    plot(x,y,'-o')
    title('u(t,x)')
    xlabel('x')
    ylabel('u(t,x)')
    axis([0,1,0,1])
    pause(0.001)
end

figure
plot(x,y,'r-o')
title('Zoom sulla regione richiesta')
axis([0.7,1,0.7,1])
xlabel('x')
ylabel('u(t,x)')

%% Le oscillazioni al bordo x=1 dipedono dal numero di Pechlet, che è minore di 1.
