clear all
close all

m=500;
h=1/(m-1);
x=linspace(0,1,m)';

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));

%condizioni ai bordi 

A(1,1:2)=[-2,2]/h^2;
A(m,m-1:m)=[0,0];
%devo finire di scrivere la condizione al bordo su u_m! La metto sul termine reazione

g=@(u) [cos(u(1:m-1));0];

%% Ora devo risolvere il problema differenziale 
y0=10*(x.^2).*(1-x);
tstar=1;
ts=4000;
k=tstar/ts;
P=phi1m(k*A);
y=y0;

for n=1:ts
	y=y+k*P*(A*y+g(y));
	%plot(x,y)
	%pause(0.1)
	%axis([0,1,-1,1])
end

yRif=y;


figure
plot(x,y)


count=0;
tsrange=10:10:100;

for ts=tsrange
count=count+1;
k=tstar/ts;
P=phi1m(k*A);
y=y0;


for n=1:ts
	y=y+k*P*(A*y+g(y));
end

err(count)=norm(y-yRif,inf);

end

figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-1),'r')
