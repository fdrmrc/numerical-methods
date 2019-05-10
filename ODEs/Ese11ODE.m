clear all
close all

A=16*[-2,1,0;1,-2,1;0,1,-2];
f=@(y) A*y + y.^2;
df=@(y) A+2*y;
y0=[sqrt(2)/2;1;sqrt(2)/2];
ts=100;
t=linspace(0,1,ts+1);
k=1/ts;
yEE=NaN(3,ts+1); %% Eulero esplicito
yEE(:,1)=y0;

for n=1:ts
  yEE(:,n+1)=yEE(:,n)+k*f(yEE(:,n));
end
figure
plot(t,yEE,'-o')
legend('Eulero esplicito')
				%Eulero implicito.

yEI=NaN(3,ts+1);
yEI(:,1)=y0;
tol =k^2;
for n=1:ts
  F=@(x) x-yEI(:,n) - k*f(x);
  J=@(x) eye(3) - k*df(x);
  x0=yEI(:,n);
  res=-J(x0)\F(x0);
  while (norm(res,Inf)>tol);
    x0+=res;
    res=-J(x0)\F(x0);
  end
  x0+=res;
  yEI(:,n+1)=x0;
end
figure
plot(t,yEI,'sk')
legend('Eulero Implicito')

				%Eulero esponenziale
yExp=NaN(3,ts+1);
yExp(:,1)=y0;

for n=1:ts
  yExp(:,n+1) = yExp(:,n)+k*phi1m(k*A)*(A*yExp(:,n)+2*yExp(:,n).^2);
end
figure
plot(t,yExp,'*')
legend('Eulero esponenziale')
