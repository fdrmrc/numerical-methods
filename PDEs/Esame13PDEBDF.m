clear all
close all

epsilon=1/100;
rho=10; %parametri del problema

h=1/100; %dato dal testo
m=150; %h=1/(m-1)
x=linspace(0,1,m)';
Pe=h/(2*epsilon)

A = toeplitz(sparse([1,2],[1,1],[-2,1]/h^2,m,1));
B = toeplitz(sparse(2,1,-1/(2*h),m,1),...
    sparse(2,1,1/(2*h),m,1));

%% condizioni al bordo (Neumann omogenee)
A(1,1:2)=[-2,2]/h^2; 
A(m,m-1:m)=[2,-2]/h^2;

B(1,1:2) = [0,0];
B(m,m-1:m) = [0,0];


%% termine reazione
g=@(u) rho*u.*(u-0.5).*(1-u);
dg=@(u) rho*(-3*u.^2+3*u-1/2);
%% ora devo rsolvere il problema differenziale y'=epsilon*A*y+B*y+reazione, con y(t0)=y0;

y0 = 10 * x.^2.*(1-x).^2 + 0.5;

%% Risoluzione problema differenziale
tsrif=1000;
ts=tsrif;
tstar=1;
k=tstar/ts;
tol=k^2/100;
y=zeros(m,ts+1); %% y=zeros(length(y0),ts+1)
t=zeros(m,ts+1);

C=epsilon*A+B;


y(:,1)=y0;
% t(1)=0;
% t(2)=t(1)+k;
y(:,2)=y0 + k*(C*y0 + g(y0));
b2=2/3;

for n=1:ts-1
    F=@(x) x - 4/3*y(:,n+1) + 1/3*y(:,n) - b2*k*(C*x+g(x));
    J=@(x) speye(m)  - b2*k*(C+spdiags(dg(x),0,m,m));
    
    yn=y(:,n);
    yn1=y(:,n+1);
    yn2=yn1;
    
    res=-J(yn2)\F(yn2);
    
    while (norm(res,inf)>tol)
    yn2=yn2+res;
    res=-J(yn2)\F(yn2);    
    end
    yn2=yn2+res;
    
    y(:,n+2)=yn2;
    
end

yrif=y;

% for j=1:ts+1
%     plot(x,y(:,j))
%     axis([0,0.2,0,1])
%       pause(0.1)
% end

count=0;
tsrange=50:50:200;

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    tol=k^2/100;
    y=zeros(m,ts+1);
    t=zeros(1,ts+1);



    
    y(:,1)=y0;
    % t(1)=0;
    % t(2)=t(1)+k;
    y(:,2)=y0 + k*(C*y0 + g(y0)); %Eulero O(k^2);
    b2=2/3;

for n=1:ts-1
    F=@(x) x - 4/3*y(:,n+1) + 1/3*y(:,n) - b2*k*(C*x+g(x));
    J=@(x) speye(m)  - b2*k*(C+spdiags(dg(x),0,m,m));
    
    yn=y(:,n);
    yn1=y(:,n+1);
    yn2=yn1;
    
    res=-J(yn2)\F(yn2);
    
    while (norm(res,inf)>tol)
    yn2=yn2+res;
    res=-J(yn2)\F(yn2);    
    end
    yn2=yn2+res;
    y(:,n+2)=yn2;
end 
   err(count)=norm(yrif(:,tsrif+1)-y(:,ts+1),inf); 
    
end


loglog(tsrange,err,'*',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-2))
