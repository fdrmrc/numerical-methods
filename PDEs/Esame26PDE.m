clear all
close all

counter=0;
mrange=2.^(4:10)+1;

for m=mrange
    counter=counter+1;
 h=(1.5*pi)/(m-1);
x=linspace(-pi,pi/2,m)';   

A = toeplitz(sparse([1,1],[1,2],[-2,1]/(h^2),1,m));
%% condizioni al bordo
A(1,1:2)=[0,0];
A(m,m-1:m)=[2,-2]/h^2;

b=@(t) 1.5*exp(t/2)*sin(x);

%% y'=Ay + b(t)
y0=sin(x);
ts=1000;
k=1/ts;
P=phi1m(k*A);
y=y0;
t=0;

for n=1:ts
    y = y+k*P*(A*y+b(t+k/2));
    t=t+k;
   % plot(x,y)
    %axis([-pi,pi/2,-2,2])
    %pause(0.01)
end
 err(counter) = norm(y-exp(t/2)*sin(x),inf);

end

loglog(mrange,err,'*',...
    mrange,err(end)*(mrange/mrange(end)).^(-2))
    
legend('errore','ordine2')
title('Error plot')
xlabel('time steps')
ylabel('error')
