clear all
close all

f=@(y) [y(2);y(3);-0.5*y(1)*y(3)];
df=@(y) [0,1,0;0,0,1;-0.5*y(3),0,-0.5*y(1)];

tstar=1;
ts=1000;
k=tstar/ts;
t=linspace(0,tstar,ts+1);
y=NaN(3,ts+1);
y0=[0,0,10]';
y(:,1)=y0;

F=@(y,yn) y -yn-k*0.5*(f(y)+f(yn));
JF=@(y) eye(3) - k*05*df(y);
tol=k^2/100;

for n=1:ts
    yn=y(:,n);
    yn1=y(:,n);  %guess iniziale
    res=-JF(yn1)\F(yn1,yn);
    
    while (norm(res,inf)>tol)
        yn1=yn1+res;
          res=-JF(yn1)\F(yn1,yn);
    end
    yn1=yn1+res;
   y(:,n+1)=yn1;
   
end

plot(t,y(1,:),'*')
disp('valore desiderato')
y(2,end)

%% uso bisezione per trovare -il valore da dare a y''(0) affinchè y'(1)=11;
a=10;
b=14;
tol=1e-5;


y=NaN(3,ts+1);
y0=[0,0,a]';
y(:,1)=y0;

for n=1:ts
    yn=y(:,n);
    yn1=y(:,n);  %guess iniziale
    res=-JF(yn1)\F(yn1,yn);
    
    while (norm(res,inf)>tol)
        yn1=yn1+res;
        res=-JF(yn1)\F(yn1,yn);
    end
    yn1=yn1+res;
    y(:,n+1)=yn1;
    
end
ya=y(2,end)-11; %% trovo il secondo estremo


y=NaN(3,ts+1);
y0=[0,0,b]';
y(:,1)=y0;
for n=1:ts
    yn=y(:,n);
    yn1=y(:,n);  %guess iniziale
    res=-JF(yn1)\F(yn1,yn);
    
    while (norm(res,inf)>tol)
        yn1=yn1+res;
        res=-JF(yn1)\F(yn1,yn);
    end
    yn1=yn1+res;
    y(:,n+1)=yn1;
    
end
yb=y(2,end)-11;


while (abs(a-b)>tol)
    
    xm=(a+b)/2;
    
    y=NaN(3,ts+1);
    y0=[0,0,xm]';
    y(:,1)=y0;    
    for n=1:ts
        yn=y(:,n);
        yn1=y(:,n);  %guess iniziale
        res=-JF(yn1)\F(yn1,yn);
        
        while (norm(res,inf)>tol)
            yn1=yn1+res;
            res=-JF(yn1)\F(yn1,yn);
        end
        yn1=yn1+res;
        y(:,n+1)=yn1;
        
    end
    ym=y(2,end)-11; %f(xn)
    
    if sign(ya)*sign(ym)<=0
        b=xm;
        yb=ym;
    else
        a=xm;
        ya=ym;
    end
    
end
disp('Il valore da dare a y''(0) è:')
xm
