clear all
close all

m=200;
a=0;
b=1;
h=(b-a)/(m-1);
x=linspace(a,b,m)';
d=1/50;


A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));

%%condizioni al bordo
A(1,1:2)=[0,0];
A(m,m-1:m)=[2,-2]/h^2;

g=@(u) [0;sin(u(2:m))];
dg=@(u) [0;cos(u(2:m))];

%% Risoluzione sistema differenziale

y0=10*x.*(1-x).^2;
yE=y0;
yER=y0;
ts=400;
tstar=1;
k=tstar/ts;

P=phi1m(k*A);
t=0;
for n=1:ts
	t+=k;
	Jn=d*A + spdiags(dg(yER),0,m,m);
	M=phi1m(k*Jn);
	yE=yE+k*P*(d*A*yE + g(yE));
	yER=yER+k*M*(d*A*yER + g(yER));
	plot(x,yE,'r-o',x,yER,'k-o')
    axis([0.95,1,0,2.5])
	  title(sprintf('Soluzione a t = %0.2f',t));
    legend('Eulero esponenziale','Eulero Rosenbrock')
    xlabel('x')
    ylabel('u(t,x)')
	pause(0.001)
end
yrifE=yE;
yrifER=yER;

count=0;
tsrange=10:10:100;

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    P=phi1m(k*A);
    yE=y0;
    yER=y0;
    
for n=1:ts
	Jn=d*A + spdiags(dg(yER),0,m,m);
	M=phi1m(k*Jn);
	yE=yE+k*P*(d*A*yE + g(yE));
	yER=yER+k*M*(d*A*yER + g(yER));
end

errE(count)=norm(yrifE - yE,inf);
errER(count)=norm(yrifER -yER,inf);

end

figure
loglog(tsrange,errE,'*',tsrange,errE(end)*(tsrange/tsrange(end)).^(-1),'r',...
    tsrange,errER,'*',tsrange,errER(end)*(tsrange/tsrange(end)).^(-2),'g')
title('Grafico logaritmico dell''errore')
xlabel('passi temporali ts')
ylabel('errori')
legend('errEulueroEsponenziale','ord1','errEuleroRosenbrock','ord2')
