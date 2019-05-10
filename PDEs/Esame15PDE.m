clear all; close all
%% DIFFERENZE FINITE + ESPONENZIALE PUNTO MEDIO
m=301;
h=pi/(m-1);
x=linspace(0,pi,m)';

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));

%% condizioni al bordo
A(1,1:2)=[0,0];
A(m,m-1:m)=[2,-2]/h^2;

b=@(t) [0;5/4*exp(t)*sin(x(2:m)/2)];

%% risolvo con esponenziale punto medio
y0=sin(x/2);

counter=0;

tsrange=[10:10:100];
for ts=tsrange
    counter=counter+1;
    k=1/ts;
    
    y=y0;
    P=phi1m(k*A);
    t=linspace(0,1,ts+1);
    
    for n=1:ts
        y= y + k*P*(A*y + b(t(n)+k/2));
    end
  
    err(counter)=norm(y- exp(1)*sin(x/2),inf);
    
end

loglog(tsrange,err,'*',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-2))
legend('error','ord1')
title('Error plot')
xlabel('timesteps')
ylabel('inf norm of error')
