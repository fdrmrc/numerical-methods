clear all
close all

m=601;
x=linspace(-1,1,m)';
h=2/(m-1);
d=1/100;c=5;

Pe=(c*h)/(2*d)<1 %verifico che il numero di Pechlet sia minore di 1

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
B = toeplitz(sparse(1,2,-1/(2*h),1,m),sparse(1,2,1/(2*h),1,m));
I=speye(m);

A(1,1:2)=[-2,2]/h^2;
B(1,2)=0;

A(m,m-1:m)=[0,0];
B(m,m-1)=0;
I(m,m)=0;
%termine reazione
%g=@(u) 0.5*[u(1:m-1);0];

y0=(x+1).^2 - 5;
tstar=0.1;
tsrif=5000;
ts=tsrif;
k=tstar/ts;
C=d*A+c*B;
y=NaN(m,ts+1);
y(:,1)=y0;
y(:,2)=(speye(m) -k*(C+0.5*I))\y(:,1);

%% coefficienti BDF 2
a0=1/3;a1=-4/3;b2=2/3;

for n=1:ts-1
    y(:,n+2)=(speye(m)-k*b2*(C+0.5*I))\(-a0*y(:,n) - a1*y(:,n+1));
end

% for j=1:ts
%     plot(x,y(:,j))
%     pause(0.001)
% end
yRef=y;
count=0;
tsrange=10:10:100;

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    
    y=NaN(m,ts+1);
y(:,1)=y0;
y(:,2)=(speye(m) -k*(C+0.5*I))\y(:,1);


for n=1:ts-1
    y(:,n+2)=(speye(m)-k*b2*(C+0.5*I))\(-a0*y(:,n) - a1*y(:,n+1));
end

err(count)=norm(y(:,ts+1)-yRef(:,tsrif+1),inf);

end

figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'r')