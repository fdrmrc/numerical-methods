clear all; close all;

%% y_n+3 - 18/11 y_n+2 +9/11y_n+1-2/11y_n=k6/11f(t_n+3,y_n+3)
%% risolvo y'=-3y, y(0)=1;
% df=@(t,y)-3;
tsrange=[10:10:100];
y0=1;
sol=@(x) exp(-3*x);
tstar=1;
counter=0;

for ts=tsrange
    counter=counter+1;    
    k=tstar/ts;
    t=linspace(0,1,ts+1);
    y=NaN(1,ts);
    y(1,1)=y0; 
    y(1,2)=sol(t(2));
    y(1,3)=sol(t(3));
    
 

for n=1:ts-2
y(1,n+3)=(18/11*y(1,n+2) - 9/11*y(1,n+1) + 2/11*y(1,n))/(1+(18*k)/11); %non serve newton visot che il problema è lineare
end

err(counter)=abs(sol(1)-y(1,ts+1));
end

plot(t,y(1,:),'sk',t,exp(-3*t),'b')

figure 
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-3),'r')
