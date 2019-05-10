clear all
close all

f=@(t,y) -0.5*y.*(1-y/3);
df=@(t,y) -0.5 + y/3;

sol=@(t) (3*exp(-0.5*t))./(2+exp(-0.5*t));
% ts=10;
% tstar=10;
% k=tstar/ts;
% t=linspace(0,tstar,ts+1);
% tol=k^3;
% 
% y=NaN(1,ts+1);
% y0=1;
% y(:,1)=y0;
% y(:,2)=sol(t(2));
% y(:,3)=sol(t(3));
% 
% %% coefficienti 
% b=6/11;
% a2=18/11; a1=9/11; a0=2/11;
% 
% F=@(tn,y,yn2,yn1,yn) y - a2*yn2 + a1*yn1 -a0*yn -k*b*f(tn+3*k,y);
% J=@(tn,y) 1 - k*b*df(tn+3*k,y);
tstar=10;
count=0;
tsrange=50:50:400;

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    t=linspace(0,tstar,ts+1);
    tol=k^3;
    
y=NaN(1,ts+1);
y0=1;
y(:,1)=y0;
y(:,2)=sol(t(2));
y(:,3)=sol(t(3));

%% coefficienti 
b=6/11;
a2=18/11; a1=9/11; a0=2/11;

F=@(tn,y,yn2,yn1,yn) y - a2*yn2 + a1*yn1 -a0*yn -k*b*f(tn+3*k,y);
J=@(tn,y) 1 - k*b*df(tn+3*k,y);
    
    
for n=1:ts-2
    yn=y(:,n);
    yn1=y(:,n+1);
    yn2=y(:,n+2);
    yn3=yn2;
    
    res=-J(t(n),yn2)\F(t(n),yn3,yn2,yn1,yn);
    while (norm(res,inf)>tol)
        yn3=yn3+res;
        res=-J(t(n),yn2)\F(t(n),yn3,yn2,yn1,yn);
    end
    yn3=yn3+res;
    
    y(:,n+3)=yn3;
    
end


err(count)=abs(y(1,ts+1)-sol(t(end)))


end

figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-3),'g')
figure
plot(t,y,'x',t,sol(t),'r')
grid on

ER=abs(y(1,end)-sol(t(end)))
