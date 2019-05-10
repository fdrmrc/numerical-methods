clear all; close all

%% Risoluzione problema differenziale 71 con RK punto medio implicito.


%% definizione parametri

a(1)=0.5;
b(1)=1;
c(1)=0.5;
g=9.81; %gravità
tsrif=4000;
ts=tsrif;
tstar=2;
y0=[pi/3,0]';
k=tstar/ts;

y=NaN(length(y0),ts+1);
y(:,1)=y0;
f=@(y) [y(2);-g*sin(y(1))];
df=@(y) [0,1;-g*cos(y(1)),0];

tol=1e-5;

for n=1:ts
    
    %% ricavo xi=w 
    w=y(:,n);
    F=@(x) x - y(:,n)-k*a(1)*f(x);
    JF=@(x) eye(2) - k*a(1)*df(x);
    
    delta=-JF(w)\F(w);
    while(norm(delta,inf)>tol)
       w=w+delta; 
       delta=-JF(w)\F(w);
    end
    
    xi(:,1)=w;
    
    y(:,n+1)=y(:,n) + k*f(xi(:,1));
    
    end
    
yrif=y;
t=linspace(0,tstar,ts+1);
plot(t,y(1,:),'*')

%% MOSTRO ORDINE DI CONVERGENZA
counter=0;

tsrange=[50:50:800];

for ts=tsrange
    
    counter=counter+1;
    k=tstar/ts;
    tol=1e-5;
    
    xi=zeros(2,1); 
    y=NaN(length(y0),ts+1);
    y(:,1)=y0;


    for n=1:ts
    
    %% ricavo xi=w 
    w=y(:,n);
    F=@(x) x - y(:,n)-k*a(1)*f(x);
    JF=@(x) eye(2) - k*a(1)*df(x);
    
    delta=-JF(w)\F(w);
    while(norm(delta,inf)>tol)
       w=w+delta; 
       delta=-JF(w)\F(w);
    end
    
    xi(:,1)=w;
    
    y(:,n+1)=y(:,n) + k*f(xi(:,1));
    
    end
    err(counter)=norm(yrif(:,tsrif+1) - y(:,ts+1),inf);

end

figure
loglog(tsrange,err,'kd',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'r',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-3),'b')
legend('errore','ordine 2','ordine 3')
