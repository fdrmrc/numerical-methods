clear all; close all;

%% Risoluzione esame 22: metodo delle linee con FD centrate + trapezi nel tempo
m=150;
h=1/(m-1);
x=linspace(0,1,m)';

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));

%% condizini al bordo
A(1,1:2)=[-2,2]/h^2;
A(m,m-1:m)=[2,-2]/h^2;

%% termine reazione 
g=@(u) u.*(u-1/2).*(1-u);
dg=@(u) ((u-1/2).*(1-u)+u.*(1-u)-u.*(u-1/2));

%% risoluzione mediante metodo dei trapezi
tsrif=1000;
ts=tsrif;
k=1/ts;
tol=k^2/100; % trapezi ha ordine 2

F=@(y,yn) y - yn - k/2*(0.02*A*(y+yn) + g(yn) + g(y));
J=@(y) eye(length(y)) - k/2*(0.02*A + spdiags(dg(y),0,m,m));

y0=10*x.^2.*(1-x).^2;
yn=y0;
for n=1:ts
    y=yn; %guess iniziale    
    delta=-J(y)\F(y,yn);    
    while (norm(delta,inf)>tol)
        y=y+delta;
        delta=-J(y)\F(y,yn);
    end
    y=y+delta;
    plot(x,y)
    axis([0,1,0,2])
    pause(0.01)
    
    yn=y;  
    
end
yRef=y;


%% convergenza 

tsrange=10:10:100;
counter=0;

for ts=tsrange
    counter=counter+1;
    k=1/ts;
    
    F=@(y,yn) y - yn - k/2*(0.02*A*(y+yn) + g(yn) + g(y));
    J=@(y) eye(length(y)) - k/2*(0.02*A + spdiags(dg(y),0,m,m));
    tol=k^2/100;
    y0=10*x.^2.*(1-x).^2;
    yn=y0;
    
    for n=1:ts
        y=yn; %guess iniziale
        delta=-J(y)\F(y,yn);
        while (norm(delta,inf)>tol)
            y=y+delta;
            delta=-J(y)\F(y,yn);
        end
        y=y+delta;
        yn=y;
        
    end
    
    err(counter)=norm(y-yRef,inf)
    
end

loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'r')
legend('errore','ordine2')
title('Error plot')
xlabel('time steps')
ylabel('error')
