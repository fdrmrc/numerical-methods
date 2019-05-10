clear all; close all

%% passaggio parametri

%a 
a(2,1)=0.5;
a(2,2)=0.5;
%b
b(1)=0.5;
b(2)=0.5;
%c
c(1)=0;
c(2)=1;

tstar=1;

y0=[1;1;1];
tsrif=3000;
ts=tsrif;
k=tstar/ts;

f=@(y) [-2*y(1)*y(2)...
    ;y(1)^2+y(3)^2-1-y(2)^2;...
    -2*(y(2)+y(1))*y(3)];df=@(y) [-2*y(2),-2*y(1),0;...
    2*y(1),-2*y(2),+2*y(3);...
    -2*y(3),-2*y(3),-2*(y(1)+y(2))];



%% inizializzo metodo
y=NaN(3,ts+1);
y(:,1) = y0(:);
xi = zeros(length(y0),2); %ho due stadi

tol=k^2/100;

for n=1:ts

    %% ricavo il primo xi
     w = y(:,n); 
    F1=@(x) x - y(:,n) - a(1,1)*k*f(x);
    JF1=@(x) eye(length(y0)) - a(1,1)*k*df(x);

  
  delta=-JF1(w)\F1(w);
  
   while (norm(delta,inf) > tol)
       w = w+delta; 
       delta = -JF1(w)\F1(w);
   end 
    w=w+delta;
    xi(:,1)=w; %trovo il primo xi1
	
    
  %% ricavo il secondo xi
  
  s= y(:,n)+a(2,1)*k*f(xi(:,1));
  
  F=@(x) x-y(:,n) - a(2,1)*k*f(xi(:,1))-a(2,2)*k*f(x);
  JF=@(x) eye(length(y0)) -a(2,2)*k*df(x);
  
  delta=-JF(s)\F(s);
  
  while(norm(delta,inf)>tol)
      s=s+delta;
      delta=-JF(s)\F(s);
  end
  s=s+delta;
  xi(:,2)=s;
  
  y(:,n+1)=y(:,n) +k*(b(1)*f(xi(:,1)) + b(2)*f(xi(:,2)));

end

yRef=y;

t=linspace(0,1,ts+1);
plot(t,y(1,:),'*')

count=0;
tsrange=[100:50:500];

for ts=tsrange
    count=count+1;
    k=tstar/ts;
    xi = zeros(length(y0),2); %ho due stadi
    y=NaN(3,ts+1);
    y(:,1)=y0;
    tol=k^2/100;

for n=1:ts

    %% ricavo il primo xi
     w = y(:,n); 
    F1=@(x) x - y(:,n) - a(1,1)*k*f(x);
    JF1=@(x) eye(length(y0)) - a(1,1)*k*df(x);

  
  delta=-JF1(w)\F1(w);
  
   while (norm(delta,inf) > tol)
       w = w+delta;
       delta = -JF1(w)\F1(w);
   end 
    w=w+delta;
    xi(:,1)=w; %trovo il primo xi1

    
  %% ricavo il secondo xi
  
  s= y(:,n)+a(2,1)*k*f(xi(:,1));
  
  F=@(x) x-y(:,n) - a(2,1)*k*f(xi(:,1))-a(2,2)*k*f(x);
  JF=@(x) eye(length(y0)) -a(2,2)*k*df(x);
  
  delta=-JF(s)\F(s);
  
  while(norm(delta,inf)>tol)
      s=s+delta;
      delta=-JF(s)\F(s);
  end
  s=s+delta;
  xi(:,2)=s;
  
  y(:,n+1)=y(:,n) +k*(b(1)*f(xi(:,1)) + b(2)*f(xi(:,2)));
end

err(count)=norm(y(:,ts+1)-yRef(:,tsrif+1),inf);
end

figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'r')
