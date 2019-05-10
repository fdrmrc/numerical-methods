clear all
close all

m=201;
x=linspace(0,1,m)';
h=1/(m-1);

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
d=1/50;

A(1,1:2)=[0,0];
A(m,m-1:m)=[2,-2]/h^2;

g=@(u) [0;cos(u(2:m))];
dg=@(u) [0;-sin(u(2:m))];

tstar=1;
tsrif=1000;
ts=tsrif;
k=tstar/ts;
tol=k^2;

F=@(y,yn) y - yn - k*(d*(A*((y+yn)/2)) + g((y+yn)/2) );
J=@(y,yn) speye(m) - k*0.5*(d*A + spdiags(dg((y+yn)/2),0,m,m) );

y=NaN(m,ts+1);
y0=x.*(x-1).^2;
y(:,1)=y0;

for n=1:ts
    yn=y(:,n);
    yn1=yn; %guess iniziale per Newton
    
    res=-J(yn1,yn)\F(yn1,yn);
    
    while (norm(res,inf)>tol)
        yn1=yn1+res;
        res=-J(yn1,yn)\F(yn1,yn);
    end
    yn1=yn1+res;
    y(:,n+1)=yn1;
    
end
t=0;
for j=1:ts+1
    t+=k;
    plot(x,y(:,j))
    title(sprintf('Soluzione a t = %0.2f',t));
    pause(0.01)
end

yRef=y;
%% ordine di convergnze
count=0;
tsrange=10:10:100;

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    
    tol=k^2;
    
    F=@(y,yn) y - yn - k*(d*(A*((y+yn)/2)) + g((y+yn)/2) );
    J=@(y,yn) speye(m) - k*0.5*(d*A + spdiags(dg((y+yn)/2),0,m,m) );
    
    y=NaN(m,ts+1);
    y0=x.*(x-1).^2;
    y(:,1)=y0;
    
    for n=1:ts
        yn=y(:,n);
        yn1=yn; %guess iniziale per Newton
        
        res=-J(yn1,yn)\F(yn1,yn);
        
        while (norm(res,inf)>tol)
            yn1=yn1+res;
            res=-J(yn1,yn)\F(yn1,yn);
        end
        yn1=yn1+res;
        y(:,n+1)=yn1;
        
    end
    
    err(count)=norm(y(:,ts+1)-yRef(:,tsrif+1),inf);
    
end

figure
loglog(tsrange,err,'g*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'r')
xlabel('ts')
ylabel('errore')
title ('errore in norma infinito')
