clear all
close all

				%Metodo utilizzato: trapezi
%N.B.: La soluzione analitica è u(x)=xsin(x)
f=@(t,y) [y(2);2*cos(x)-y(1)];
y0=[0;0];
A=[0,1;-1,0];
tsrange=[20:10:300];
count=0;

for ts=tsrange
  count++;
  k=1/ts;
  y=NaN(2,ts+1);
  y(:,1)=y0;

  t=linspace(0,1,ts+1);
  b=@(t)[0;2*cos(t)];
  
  for n=1:ts
    y(:,n+1)=(eye(2)-0.5*k*A)\(y(:,n)+k*0.5*(A*y(:,n)+b(t(n))+b(t(n+1))) );
  endfor
  
  err(count)=norm(y(1,:)-t.*sin(t),Inf); 
endfor


loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'g')
title('Error plot')
xlabel('time steps')
ylabel('inf norm of error')
  
