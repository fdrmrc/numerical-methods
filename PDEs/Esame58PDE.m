clear all
close all

m=150;
h=1/(m-1);
x=linspace(0,1,m)';

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
B = toeplitz(sparse(1,2,-1/(2*h),1,m), sparse(1,2,1/(2*h),1,m));

% condizioni al bordo

A(1,1:2)=[0,0];
B(1,2)=0;

A(m,m-1:m)=[2,-2]/h^2;
B(m,m-1:m)=[0,0];

% risoluzione problema mediante eulero-Rosenbrock esponenziale
y0=x.*(1-x).^2;
tstar=1;
ts=1000;
k=tstar/ts;

y=y0;
for n=1:ts

	Jn=A+2*(B^2)*diag(y);
	P=phi1m(k*Jn);
	y=y + k*P*(A*y + (B^2)*y.^2);
	%plot(x,y)
	%pause(0.1)
	
end

%figure
%plot(x,y)
yrif=y;

count=0;
tsrange=10:10:100;

for ts=tsrange

	count=count+1;	
	k=tstar/ts;
	y=y0;

for n=1:ts

	Jn=A+2*(B^2)*diag(y);
	P=phi1m(k*Jn);
	y=y + k*P*(A*y + (B^2)*y.^2);
end

err(count)=norm(y-yrif,inf);
end

figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'g')
title('error plot')
xlabel('timesteps')
ylabel('error w.r.t. sol riferimento')
legend('error','ord2')
