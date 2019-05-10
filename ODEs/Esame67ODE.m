clear all
close all

f=@(t,y) [y(2);y(3);-100*y(3)-t*y(2)-t^2*y(1)*y(2)];
df=@(t,y)[0,1,0;0,0,1;-y(2)*t^2,-t-y(1)*t^2,-100];

tstar=10;
tsrif=1000;
ts=tsrif;
k=tstar/ts;
t=linspace(0,tstar,ts+1);
y0=[1,2,0]';
y=NaN(length(y0),ts+1);
y(:,1)=y0;
tol=k^2/100;
%Scelta: trapezi. Infatti è A-Stabile!
F=@(tn,y,yn) y-yn-k*0.5*(f(tn,yn) + f(tn+k,y));
J=@(tn,y) eye(3) - k*0.5*df(tn+k,y);

for n=1:ts
    yn=y(:,n);
    yn1=yn; %guess iniziale
    
    res=-J(t(n),yn1)\F(t(n),yn1,yn); 
    while (norm(res,inf)>tol)
        yn1=yn1+res;
        res=-J(t(n),yn1)\F(t(n),yn1,yn); 
    end
    yn1=yn1+res; %iterazione extra, ho già trovato un residuo in più.
    
    y(:,n+1)=yn1;
    
end

plot(t,y(1,:),'*')
yRef=y;
tsrange=10:10:100;
count=0;

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    t=linspace(0,tstar,ts+1);
    y=NaN(3,ts+1);
    y(:,1)=y0;
    
    
    F=@(tn,y,yn) y-yn-k*0.5*(f(tn,yn) + f(tn+k,y));
J=@(tn,y) eye(3) - k*0.5*df(tn+k,y);

for n=1:ts
    yn=y(:,n);
    yn1=yn; %guess iniziale
    
    res=-J(t(n),yn1)\F(t(n),yn1,yn); 
    while (norm(res,inf)>tol)
        yn1=yn1+res;
        res=-J(t(n),yn1)\F(t(n),yn1,yn); 
    end
    yn1=yn1+res; %iterazione extra, ho già trovato un residuo in più.
    
    y(:,n+1)=yn1;
    
end

err(count)=norm(y(:,ts+1)-yRef(:,tsrif+1),inf);

end

figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'r')

%% Verifico se si tratta di un problema stiff

cStiff=0; %count stiff
t=linspace(0,tstar,tsrif+1);

for j=1:tsrif
    v=eig(df(t(j),yRef(:,j)));
    if real( (v(1))<0 && real(v(2)<0) ) && (abs(v(1)-v(2)) > 10)
        cStiff=cStiff+1;
    end
end

cStiff
