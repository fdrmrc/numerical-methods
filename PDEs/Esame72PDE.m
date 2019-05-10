clear all
close all

m=501;
a=0;
b=1;
c=5/2;
d=1/100;

h=(b-a)/(m-1);
x=linspace(a,b,m)';

Pe=(c*h)/(2*d)
Pe<1


A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
B = toeplitz(sparse(1,2,-1/(2*h),1,m),sparse(1,2,1/(2*h),1,m));

%% condizioni al bordo
A(1,1:2)=[-2,2]/h^2;
%A(1,1:2)=[0,0];
B(1,1:2)=[0,0];

A(m,m-1:m)=[0,0];
B(m,m-1:m) = [0,0];
I=ones(m,1);
I(m,1)=0;

b=@(t) I*cos(2*t);  % in questo modo l'ultima componente è nulla


%% risoluzione problema differenziale
y0=(x.^2).*(1-x);
tstar=1/2;
tsrif=2000;
ts=tsrif;
k=tstar/ts;
C=d*A+c*B;
P=phi1m(k*C);
t=0;
y=y0;


for n=1:ts
    y=y + k*P*(C*y + b(t+k/2));
    t=t+k;
   	plot(x,y)
    title(sprintf('Soluzione a t = %0.2f',t));
    axis([-0.1,0.1,0,0.5]) %zoommo sul bordo a sx
    pause(0.01)
    
end
yrif=y;

counter=0;
tsrange=10:10:100;

for ts=tsrange
    counter=counter+1;
    k=tstar/ts;    
    C=d*A+c*B;
    P=phi1m(k*C);
    t=0;
    y=y0;    
    
    for n=1:ts
        y=y + k*P*(C*y + b(t+k/2));
        t=t+k;
    end
    
    err(counter)=norm(y- yrif,inf);
    
end

loglog(tsrange,err,'*',tsrange,(err(end))*(tsrange/tsrange(end)).^-2,'r')
legend('error','ordine2')