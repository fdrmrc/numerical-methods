clear all
close all

count=0;
mrange=200:50:900;


for m=mrange 
    count=count+1;    
    x=linspace(0,1,m)';
h=1/(m-1);

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
B = toeplitz(sparse(1,2,-1/(2*h),1,m), sparse(1,2,1/(2*h),1,m));

c=150;
d=4/(pi^2);
Pe=(c*h)/(2*d);
Pe<1

%phi=@(x) x;
phiS=@(z) z-1+ (2*z)./(exp(2*z)-1) ;
%d=d*(1+phi(Pe));
d=d*(1+phiS(Pe));

%condizioni al bordo
A(1,1:2)=[0,0];
B(1,2)=0;
A(m,m-1:m)=[2,-2]/h^2;
B(m,m-1)=0;
%termine sorgente
b=@(t) c*(pi/2)*exp(-t)*[0;cos((pi/2)*x(2:m))];

ts=1000;
tstar=1;
k=tstar/ts;
y0=sin((pi/2)*x);
t=linspace(0,tstar,ts+1);
C=d*A-c*B;
yn=y0;

for n=1:ts
    y=yn;
    y=(speye(m)-k*0.5*C)\( (speye(m)+k*0.5*C)*y + k*0.5*( b(t(n)) + b(t(n+1)) ) );
    yn=y;
end

sol=exp(-1)*sin((pi/2)*x);
err(count)=norm(y-sol,inf);

end

loglog(mrange,err,'*',mrange,err(end)*(mrange/mrange(end)).^(-1),'g',...
    mrange,err(end)*(mrange/mrange(end)).^(-2),'k')
title('Spatial error plot')
xlabel('m')
ylabel('Error')
legend('error','ord1','ord2')
