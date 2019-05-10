clear all; close all

m=200;
d=1/100;
x = linspace(0,1,m)';
h = 1/(m-1);
A = toeplitz(sparse([1,1],[1,2],... 
    [-2,1]/h^2,1,m)); %discretizzazione spaziale
g=@(u) 5*u.*(u-1).*(1-u);
dg= @(u) -2.5*(6*u.^2-6*u+1);

%% condizioni al bordo
A(1,1:2)=[0,0];
A(m,m-1:m)=[0,0];

% Devo azzerare la prima riga del termine "reazione".
ts=1000;
tstar=1;
k=tstar/ts;
tol=k^2;
u0=4*x.*(1-x);
yn=u0; %valore iniziale
t=0;

for n=1:ts
    y=yn; %guess iniziale
    
    F=@(y,yn) y -yn -k*d*A*y - k*[0;g(y(2:m-1));0];
    JF=@(y) eye(length(y)) - k*d*A - k*diag([0;dg(y(2:m-1));0]);
    
    delta=-JF(y)\F(y,yn);	
    while(norm(delta,inf)>tol)
        y=y+delta;
        delta=-JF(y)\F(y,yn);
    end
    y=y+delta;
    t+=k;
    
        plot(x,y,'-o')
        title(sprintf('Soluzione a t = %0.2f',t));
	xlabel('x')
	ylabel('u(t,x)')
        axis([0,1,-0.3,1.3])
        pause(0.0001)
    yn=y; %gli passo quello appena trovato (atrimenti riprende yn=u0)
end
yrif=y;
count=0;
tsrange=10:10:50;

for ts=tsrange
    count+=1;
    k=tstar/ts;
    tol=k^2/100;
    yn=u0;
    
    for n=1:ts
        y=yn; %guess iniziale
        
        F=@(y,yn) y -yn -k*d*A*y - k*[0;g(y(2:m-1));0];
        JF=@(y) eye(length(y)) - k*d*A - k*diag([0;dg(y(2:m-1));0]);
        
        delta=-JF(y)\F(y,yn);
        while(norm(delta,inf)>tol)
            y=y+delta;
            delta=-JF(y)\F(y,yn);
        end
        y=y+delta;
        
        yn=y; 
    end
    err(count)=norm(y-yrif,inf);
end

figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-1),'g')
title('Error plot')
legend('inf norm errore','ordine 1') 
