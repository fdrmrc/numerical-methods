clear all
close all

m=501;
tstar=1;
a=pi/4;
b=pi/2;
x=linspace(a,b,m)';
h=(b-a)/(m-1);

A = toeplitz(sparse([1,1],[1,2],[-2,1]/(h^2),1,m));
d=1/4;
A(1,1:2)=[-2,2]/h^2;
A(m,m-1:m)=[0,0];

tsrange=10:10:100;
count=0;
for ts=tsrange
    count=count+1;
    k=tstar/ts;
    tol=k^2;
    y0=sin(2*x);
    t=linspace(0,tstar,ts+1);
    y=NaN(m,ts+1);
    y(:,1)=y0;
    
    F=@(y,yn,t) y - yn - k*(d*(A*y) - [y(1:m-1);0].^2 + exp(-2*t)*(sin(2*x)).^2);
    J= @(y) speye(m) - k*(d*A - diag(2*[y(1:m-1);0]));
   
    
    for n=1:ts
        yn=y(:,n);
        yn1=yn;
        res=-J(yn1)\F(yn1,yn,t(n+1));
        while (norm(res,inf)>tol)
            yn1=yn1+res;
            res=-J(yn1)\F(yn1,yn,t(n+1));
        end
        yn1=yn1+res;
        y(:,n+1)=yn1;
    end
    
    err(count)=norm(y(:,ts+1)-exp(-1)*sin(2*x),inf);
    
end

loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-1),'g')