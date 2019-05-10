clear all
close all

m=501;
x=linspace(0,pi,m)';
h=pi/(m-1);

 A= toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
 I=speye(m);
 
 %condizioni al bordo
 A(1,1:2)=[0,0];
 A(m,m-1:m)=[2,-2]/h^2;
 I(1,1)=0;
 b=@(t) x*sin(t);
 
 %% Risoluzione sistema differenziale
 y0=sin(x/2);
 tsrif=1000;
 ts=tsrif;
 tstar=1;
 k=tstar/ts;
 y=NaN(length(y0),ts+1);
 y(:,1)=y0;
 t=linspace(0,tstar,ts+1);
  C=A-3*I;
  P=phi1m(k*C);
  
 for n=1:ts
     
     y(:,n+1)=y(:,n)+k*P*(C*y(:,n) + b(t(n)+k/2) );
 end
 yRef=y; %% SOLUZIONE DI RIFERIMENTO
 t=0;
 for j=1:ts
 t+=k;
     plot(x,y(:,j))
     	  title(sprintf('Soluzione a t = %0.2f',t));
     axis([0,pi,0,0.9])
     pause(0.001)
 end
 
 count=0;
 tsrange=10:10:100;
 
 for ts=tsrange
     count=count+1;
     k=tstar/ts;
        
     y=NaN(length(y0),ts+1);
    y(:,1)=y0;
    t=linspace(0,tstar,ts+1);
    C=A-3*I;
    P=phi1m(k*C);
  
 for n=1:ts
     
     y(:,n+1)=y(:,n)+k*P*(C*y(:,n) + b(t(n)+k/2) );
 end
 
 err(count)=norm(y(:,ts+1)-yRef(:,tsrif+1),inf);
 
 end
 
 figure
 loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2))
