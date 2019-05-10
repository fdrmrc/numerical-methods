clear all; close all

m=150;
h=1/(m-1);
x=linspace(0,1,m)';

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));

%% condizioni al bordo
A(1,1:2)=[0,0];
A(m,m-1:m)=[2,-2]/h^2;


%% risoluzione problema differenziale (con metodo dei trapezi)
y0=10 * x.*(1-x).^2;
ts=1000;
k=1/ts;

F=@(y,yn) y -yn - k/2*(0.01*A*yn + [0;sin(yn(2:m))]) - k/2*(0.01*A*y + [0;sin(y(2:m))]);
JF=@(y) eye(length(y)) - k/2*(0.01*A + spdiags([0;cos(y(2:m))],0,m,m));
yn=y0;
tol=k^2/100;

for n=1:ts
    y=yn;
    delta=-JF(y)\F(y,yn);
    while(norm(delta,inf)>tol)
        y=y+delta;
        delta=-JF(y)\F(y,yn);
    end
    y=y+delta;
    plot(x,y)
  %  axis([0,1,-1,5])
    pause(0.001)
    yn=y;
end

yRif=y;
plot(x,y)
counter=0;
tsrange=[10:10:100];
y0=10 * x.*(1-x).^2;

for ts=tsrange
    counter=counter+1;
 
    k=1/ts;
    F=@(y,yn) y -yn - k/2*(0.01*A*yn + [0;sin(yn(2:m))]) - k/2*(0.01*A*y + [0;sin(y(2:m))]);
    JF=@(y) eye(length(y)) - k/2*(0.01*A + spdiags([0;cos(y(2:m))],0,m,m));    yn=y0;
    tol=k^2/100;
    
    for n=1:ts
        y=yn;
        delta=-JF(y)\F(y,yn);
        while(norm(delta,inf)>tol)
            y=y+delta;
            delta=-JF(y)\F(y,yn);
        end
        y=y+delta;
        yn=y;
    end
    
    err(counter)=norm(y-yRif, inf);
    
end
loglog(tsrange,err,'*',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-1),'r',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'b')

legend('errore','ordine 1','ordine 2')
