clear all
close all

mrange=[10:10:100];
count=0;

for m=mrange+1
    count=count+1;
    x=linspace(1,2,m)';
    h=1/(m-1);
    D2 = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
    D1 = toeplitz(sparse(1,2,-1/(2*h),1,m), sparse(1,2,1/(2*h),1,m));

    
    
    %condizioni al bordo
    b=zeros(m,1);
    
    D2(1,1:2)=[1,0];
    D1(1,1:2)=[0,0];
    
    D2(m,m-1:m)=[0,1];
    D1(m,m-1:m)=[0,0];
    
    b(m,1)=-1/2;
    
    F=@(u) D2*u - D1*((1-u).^2) + b;
    J=@(u) D2 + 2*D1 - 2*spdiags([[-u(2:m-1)/h;0;0],[0;(u(3:m)-u(1:m-2))/h;0],[0;0;u(2:m-1)/h]],[-1,0,1],m,m);
    
    u=ones(m,1);
    
      res=-J(u)\F(u);
    
    tol=h^2/100;
    while (norm(res,inf)>tol)
        u=u+res;
        res=-J(u)\F(u);
    end
    u=u+res;
    sol=1-1./x;
    err(count)=norm(u-sol,inf);
    
end

figure
plot(x,u,'*',x,1-1./x,'g')
xlabel('x')
ylabel('y')
title('Soluzione')
legend('soluzione numerica','sol analitica')

figure
loglog(mrange,err,'*',mrange,(err(end))*(mrange/mrange(end)).^-2)
xlabel('m')
ylabel('errore')    
