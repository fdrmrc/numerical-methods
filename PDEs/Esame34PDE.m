clear all
close all

m=51;
h=1/(m-1);
x=linspace(0,1,m)';

A = toeplitz(sparse([1,1],[1,2],[-2,1]/(h^2),1,m)); 
B = toeplitz(sparse(1,2,-1/(2*h),1,m), sparse(1,2,1/(2*h),1,m));

%condizioni al bordo
A(1,1:2)=[-2,2]/h^2;
B(1,2)=0;

A(m,m-1:m)=[0,0];
B(m,m-1)=0;

d=4;c=2;
C=d*A-c*B;

y0=x.^2;
tstar=1;
tsrif=1e5;
ts=tsrif;
k=tstar/ts
y=y0;

for n=1:ts
    y=y+k*(C*y);
%     plot(x,y)
%     pause(0.001)
end
yRef=y;
tsrange=20000:10:20100;
count=0;

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    
    y=y0;

for n=1:ts
    y=y+k*(C*y);
end

err(count)=norm(y-yRef,inf);

end

loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-1),'b',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'r')

%% La restrizione sul passo temporale è dovuta al fatto che eulero esplicito non è A-stabile, dunque serve una restrizione. Per calcolare quanto deve essere, basta osservare a cosa si riduce il sistema "semidiscretizzato" (cioè dopo aver applicato differenze finite agli operatori spaziali).
