clear all; close all

c=4;d=0.1;p=50;

m=101;
h=1/100;
x=linspace(0,1,m)';

A = toeplitz(sparse([1,1],[1,2],[-2,1]/(h^2),1,m));
B = toeplitz(sparse(1,2,-1/(2*h),1,m),...
    sparse(1,2,1/(2*h),1,m));

% condizioni al bordo
A(1,1:2)=[-2,2]/h^2;
B(1,2)=0;

A(m,m-1:m)=[0,0];
B(m,m-1)=0;

g=@(u) p*u.^2.*(1-u);
dg=@(u) p*(2*u - 3*u.^2);

y0=-x.^2 +1;
ts=1000;
tstar=0.1;
k=tstar/ts;

F=@(y,yn) y- yn- k/2*(d*A*yn - c*B*yn + [g(yn(1:m-1));0] + d*A*y - c*B*y + [g(y(1:m-1));0]);
J=@(y) eye(length(y)) - k/2*(d*A-c*B+spdiags([dg(y(1:m-1));0],0,m,m));
yn=y0;
tol=k^2/100;

for n=1:ts
    y=yn;
    delta=-J(y)\F(y,yn);
    
    while (norm(delta,inf)>tol)
        y=y+delta;
        delta=-J(y)\F(y,yn);
    end
    y=y+delta;
    
    yn=y;
end
yRef=y;
%plot(x,y)

tsrange=10:10:100;
counter=0;

for ts=tsrange
    counter=counter+1;
    k=tstar/ts;
    
    F=@(y,yn) y- yn- k/2*(d*A*yn - c*B*yn + [g(yn(1:m-1));0] + d*A*y - c*B*y + [g(y(1:m-1));0]);
    J=@(y) eye(length(y)) - k/2*(d*A-c*B+spdiags([dg(y(1:m-1));0],0,m,m));
    yn=y0;
    tol=k^2/100;
    
    for n=1:ts
        y=yn;
        delta=-J(y)\F(y,yn);
        
        while (norm(delta,inf)>tol)
            y=y+delta;
            delta=-J(y)\F(y,yn);
        end
        y=y+delta;
        
        yn=y;
    end
    
    err(counter)=norm(y-yRef,inf);
    
end
    
loglog(tsrange,err,'*',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-2))
