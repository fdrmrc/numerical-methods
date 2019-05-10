clear all
close all

tstar=1;
tsrif=1000;
ts=tsrif;
k=tstar/ts;

f=@(t,y) -10*(y-cos(t));
df=@(t,y) -10; % IL PROBLEMA E' LINEARE=> NON SERVE NEWTON
y0=0;
y=NaN(1,ts+1);
t=linspace(0,tstar,ts+1);
y(1,1)=y0;
%y(1,2)=y0+k*(f(t(2),y0)); %O(k^2) => Poi metti trapezi
%y(1,2)= 10/101*(-10*exp(-10*t(2)) + 10*cos(t(2)) + sin(t(2))); %% In
%questo caso SPECIALE conosco la soluzione. 
%Altrimenti, dovrei usare un metodo che LOCALMENTE ha ordine 3=> trapezi.
%Il problema resta lineare perciò devono solamente risolvere un sistema
%LINEARE.
y(1,2)=(1+5*k)\( (1-5*k)*y(1,1)+5*k*cos(t(1)) + 5*k*cos(t(2)) );

for n=1:ts-1
    y(1,n+2)=(1+k*(50/12))\...
        ((1-k*(80/12))*y(1,n+1)+10*(k/12)*y(1,n) + k*(50/12)*cos(t(n+2)) + k*(20/3)*cos(t(n+1)) - (5/6)*k*cos(t(n)));
end

yRef=y;
sol=@(x) 10/101*(-10*exp(-10*t) + 10*cos(t) + sin(t));

plot(t,y(1,:),'*',t,sol(t),'k')

count=0;
tsrange=[50:50:200];

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    t=linspace(0,tstar,ts+1);
    y0=0;
    
y=NaN(1,ts+1);
t=linspace(0,tstar,ts+1);
y(1,1)=y0;
%y(1,2)= 10/101*(-10*exp(-10*t(2)) + 10*cos(t(2)) + sin(t(2)));
y(1,2)=(1+5*k)\( (1-5*k)*y(1,1) +5*k*cos(t(1)) + 5*k*cos(t(2)) );
for n=1:ts-1
    y(1,n+2)=(1+k*(50/12))\...
        ((1-k*(80/12))*y(1,n+1)+10*(k/12)*y(1,n) + k*(50/12)*cos(t(n+2)) + k*(20/3)*cos(t(n+1)) - (5/6)*k*cos(t(n)));
end
err(count)=norm(y(1,ts+1)-yRef(1,tsrif+1),inf);

end

figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'r',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-3),'g')
legend('Errore','Ordine 2','Ordine 3')
