clear all
close all

m=601;
h=(pi/2)/(m-1);
x=linspace(0,pi/2,m)';

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
A(1,1:2)=[-2,2]/h^2;
A(m,m-1:m)=[0,0]; 

%termine sorgente
b=@(t) (cos(x)/2)*(cos(t/2) - sin(t/2));
%non va azzerata pure ultima riga perché il nodo xm=pi/2 la azzera
%automaticmaente

%% Risoluzione sistema differenziale
y0=cos(x);
tstar=1;
d=0.5;
t=0;

tsrange=10:10:100;
count=0;
for ts=tsrange
    count=count+1;
    k=tstar/ts;
    y=NaN(m,ts+1);
    
    y(:,1)=y0;
    t=0;

    
    P=phi1m(k*d*A);
for n=1:ts
    y(:,n+1)=y(:,n) + k*P*(d*A*y(:,n) +b(t+k/2));
    t=t+k;
end

err(count)=norm(y(:,ts+1)-cos(x)*cos(1/2),inf);

end

loglog(tsrange,err,'*',tsrange,(err(end))*(tsrange/tsrange(end)).^(-2) )

% 
% for j=1:ts
%     plot(x,y(:,j))
%     pause(0.001)
% end



% yRef=y;
% %% ordine di convergenza temporale
% count=0;
% tsrange=10:10:100;
% 
% for ts=tsrange
%     count=count+1;
%     k=tstar/ts;
%     
%     P=phi1m(k*d*A);
%     t=0;
%     for n=1:ts
%         y(:,n+1)=y(:,n) + k*P*(d*A*y(:,n) +b(t+k/2));
%         t=t+k;
%     end
%     err(count)=norm(y(:,ts+1)-yRef(:,tsrif+1),inf)
%     
% end
% figure
% loglog(tsrange,err,'*',tsrange,(err(end))*(tsrange/tsrange(end)).^(-2) )