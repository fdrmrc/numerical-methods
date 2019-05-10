clear all; close all

%% Risoluzione esame PDE numero 20: metodo delle linee: FD centrate + Eulero implicito
m=501;
h=1/(m-1);
x=linspace(0,1,m)';

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));

% condizioni al bordo
A(1,1:2)=[-2,2]/h^2;
A(m,m-1:m)=[0,0];

%devo ancora finire di mettere le condizioni al bordo, poiché ho un termine
%di reazione g
g=@(u) [cos(u(1:m-1));0];

% mi sono ricondotto al sistema y'=0.5*A*y + g

%% risoluzione sistema differenziale usando eulero implicito.
tsrif=1000;
ts=tsrif;
k=1/ts;
y0=10*x.^2.*(1-x) + 1;

F=@(y,yn) y- yn - k*(0.02*A*y + [cos(y(1:m-1));0]);
J=@(y) eye(length(y)) - k*0.02*A - k*spdiags([-sin(y(1:m-1));0],0,m,m);
yn=y0;
tol=k^2/100;

for n=1:ts
    y=yn;
    delta=-J(y)\F(y,yn);
    while (norm(delta,inf)>tol)
        y=y+delta;
        delta=-J(y)\F(y,yn);
    end
    y=y+delta;
    plot(x,y)
    pause(0.01)
    yn=y;
end

yRif=y;
max(yRif)
%% convergenza
counter=0;
tsrange=10:10:100;

for ts=tsrange
    counter=counter+1;
    k=1/ts;
    
    F=@(y,yn) y- yn - k*(0.02*A*y + [cos(y(1:m-1));0]);
    J=@(y) eye(length(y)) - k*0.02*A -k*spdiags([-sin(y(1:m-1));0],0,m,m);
    yn=y0;
    tol=k^2/100;
    
    for n=1:ts
        y=yn;
        delta=-J(y)\F(y,yn);
        while (norm(delta,inf)>tol)
            y=y+delta;
            delta=-J(y)\F(y,yn);
        end
        y=y+delta;
        yn=y;
    end
    err(counter)=norm(y - yRif, inf);    
end

loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-1))
title('Plot dell''errore')
legend('errore','ordine1')
