clear all
close all

m=301;
h=1/(m-1);
x=linspace(0,1,m)';
d=1/25;
A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));

A(1,1:2)=[-2,2]/h^2;
A(m,m-1:m)=[0,0];
I=speye(m);
I(m,m)=0;

tstar=1;
ts=3000;
k=tstar/ts;

y0=(x.^2).*(1-x);
yn=y0;
tol=k^2/100;

for n=1:ts
    
    y=yn;
    F=@(y,yn) y - yn - k*(d*A*y + ...
        I*spdiags(x,0,m,m)*(1./(1+y.^2)));
    J=@(y) eye(length(y)) - k*d*A - k*I*spdiags(x,0,m,m)*spdiags((-2*y)./((1+y.^2).^2),0,m,m);
    
    res=-J(y)\F(y,yn);
    while (norm(res,inf)>tol)
        y=y+res;
        res=-J(y)\F(y,yn);
    end
    y=y+res;
    
    yn=y;
    
    
end


yrif=y;

count=0;
tsrange=10:10:100;

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    yn=y0;
    tol=k^2/100;
    
    for n=1:ts
        
        y=yn;
        F=@(y,yn) y - yn - k*(d*A*y + ...
            I*spdiags(x,0,m,m)*(1./(1+y.^2)));
        J=@(y) eye(length(y)) - k*d*A - k*I*spdiags(x,0,m,m)*spdiags((-2*y)./((1+y.^2).^2),0,m,m);
        
        res=-J(y)\F(y,yn);
        while (norm(res,inf)>tol)
            y=y+res;
            res=-J(y)\F(y,yn);
        end
        y=y+res;
        
        yn=y;
        
        
    end
    err(count)=norm(y-yrif,inf);
    
end

loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-1),'r')
title('Error plot')
legend('error','ord1')