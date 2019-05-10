clear all; close all

%% Risoluzione problema differenziale numero 75, con metodo dei trapezi

tstar=1;
tsrif=4000;
ts=tsrif;
k=tstar/ts;

y0=[1,0,1]';
theta=0.5;

y=NaN(length(y0),ts+1);
y(:,1)=y0;

%% definizione funzione
f=@(t,y) [y(1)+y(2)-cos(t);y(3);y(3)+(y(2).^2).*y(1)];
df=@(t,y) [1,1,0;0,0,1;1+y(2).^2,2*y(1).*y(2),0];

F = @(tn,y,yn) y-yn-k/2*(f(tn,yn)+f(tn+k,y));
JF = @(tn,y) eye(length(y0))-k/2*df(tn+k,y);
t=linspace(0,tstar,ts+1);
tol=k^2/100;

for n=1:ts
    
    y(:,n+1)=y(:,n);
    
    delta=-JF(t(n),y(:,n+1))\F(t(n),y(:,n+1),y(:,n));
    
    while(norm(delta,inf)>tol)
        y(:,n+1)=y(:,n+1)+delta;
        delta=-JF(t(n),y(:,n+1))\F(t(n),y(:,n+1),y(:,n));
    end
    y(:,n+1)=y(:,n+1)+delta;
    
end
yrif=y;
 
 plot(t,y)
 counter=0;
 
 tsrange=[50:50:500];
 for ts=tsrange
     counter=counter+1;
     k=tstar/ts;
     y=NaN(length(y0),ts+1);
     y(:,1)=y0;
     tol=k^2/100;
     t=linspace(0,tstar,ts+1);
     
     F = @(tn,y,yn) y-yn-k/2*(f(tn,yn)+f(tn+k,y));
     JF = @(tn,y) eye(length(y0))-k/2*df(tn+k,y); 
     for n=1:ts
         
         y(:,n+1)=y(:,n);
         
         delta=-JF(t(n),y(:,n+1))\F(t(n),y(:,n+1),y(:,n));
         
         while(norm(delta,inf)>tol)
             y(:,n+1)=y(:,n+1)+delta;
             delta=-JF(t(n),y(:,n+1))\F(t(n),y(:,n+1),y(:,n));
         end
         y(:,n+1)=y(:,n+1)+delta;
     end
     
     err(counter)=norm(yrif(:,tsrif+1) - y(:,ts+1),inf);
     
 end
 
 
 figure
 loglog(tsrange,err,'k*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'r',...
     tsrange,err(end)*(tsrange/tsrange(end)).^(-3),'b')
 legend('errore','ordine 2','ordine 3')
 
% %%Check
F(0,y(:,2),y0)
