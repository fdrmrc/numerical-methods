clear all
close all

m=600;
h=1/(m-1);
x=linspace(0,1,m)';

c=+1;d=0.001;
Pe=(abs(c)*h)/(2*d)

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
B = toeplitz(sparse(1,2,-1/(2*h),1,m), sparse(1,2,1/(2*h),1,m));

%condizioni al bordo 
A(1,1:2)=[-2,2]/h^2;
B(1,1:2)=[0,0]; 	%Neumann

A(m,m-1:m)=[0,0];
B(m,m-1:m)=[0,0];
C=d*A+c*B;

V=ones(m,1);
V(m,1)=0;
b=@(t) - V.*log(1+2*t);

y0=(x.^2).*(1-x);
tsrif=1000;
ts=tsrif;
tstar=1;
k=tstar/ts;
P=phi1m(k*C);
t=0;
y=y0;

for n=1:ts	
	y=y + k*P*(C*y + b(t + k/2));
	t=t+k;
	plot(x,y)
        title(sprintf('Soluzione a t = %0.2f',t));
     axis([0,1,-1,1])
	pause(0.01)		
end

yrif=y;
count=0;
tsrange=10:10:200;

for ts=tsrange
	count=count+1;
	k=tstar/ts;
    
    P=phi1m(k*C);
    t=0;
    y=y0;

for n=1:ts	
	y=y + k*P*(C*y + b(t + k/2));
	t=t+k;
end


err(count)=norm(y-yrif,inf);

end

figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'r')
title('Grafico dell''errore')
xlabel('timesteps')
ylabel('errore')
legend('errore','ord2')
