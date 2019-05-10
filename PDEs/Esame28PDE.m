clear all; close all

%% Risoluzione esame 28 con FD + Eulero implicito

m=101;
h=1/(m-1);
x=linspace(0,1,m)';

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));

%% condizioni al bordo
A(1,1:2)=[0,0];
A(m,m-1:m)=[0,0]; 

counter=0;
tsrange=10:10:100;
for ts=tsrange
    counter=counter+1;
    k=1/ts;
    y0=x.^2-x;
y=y0;
c=(x-x.^2)/2;

for n=1:ts
   y=(eye(length(y0))-k*diag(c)*A)\y; 
end
err(counter)=norm(y-exp(-1)*(x.^2-x),inf);
%err(counter)=norm(y-yRef,inf);

    
end
figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-1),'r')

%% ordine spaziale è costante perché la soluzione al tempo 1 è un polinomio di grado 2 e ha derivata quarta nulla
count=0;

