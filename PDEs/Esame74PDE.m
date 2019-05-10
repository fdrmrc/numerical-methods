clear all
close all

m=600;
h=1/(m-1);
x=linspace(0,1,m)';

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
B = toeplitz(sparse(1,2,-1/(2*h),1,m), sparse(1,2,1/(2*h),1,m));
d=1/100;
c=5;

%%%PECHLET!
Pe=(c*h)/(2*d)


A(1,1:2)=[0,0];
B(1,2)=0;
A(m,m-1:m)=[2,-2]/h^2;
B(m,m-1)=0;

%% termine noto
V=ones(m,1);
V(1,1)=0;
b=@(t) V.*sin(t);

% risoluzione con eulero implicito

y0=x.*(1-x).^2;
tstar=0.1;
ts=2000;
k=tstar/ts;
C=d*A-c*B;
t=linspace(0,tstar,ts+1);
y=y0;
for n=1:ts
    y=(speye(m)-k*C)\(y+k*b(t(n+1)));
%     plot(x,y)
%     pause(0.1)
end
yrif=y;


tsrange=50:50:400;
count=0;

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    t=linspace(0,tstar,ts+1);
    y=y0;
    %[L,U,P]=lu(speye(m)-k*C);
    
for n=1:ts
   y=(speye(m)-k*C)\(y+k*b(t(n+1)));
   % y=U\(L\(P*(k*b(t(n+1)) + y)));
end

err(count)=norm(y-yrif,inf);

end

figure
loglog(tsrange,err,'*',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-1),'r',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'g')
legend('error','ord1','ord2')