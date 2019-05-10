clear all
close all
tic

f=@(t,y) [y(2)*y(3)*sin(t)-y(1)*y(2)*y(3);...
	  -y(1)*y(3)*sin(t)+(1/20)*y(1)*y(3);...
	  y(1)^2*y(2) - (1/20)*y(1)*y(2)];

df=@(t,y) [y(2)*y(3),y(3)*sin(t)-y(1)*y(3),y(2)*sin(t)-y(1)*y(2);...
	   -y(3)*sin(t)+(1/20)*y(3),0,-y(1)*sin(t)+(1/20)*y(1);...
	   2*y(1)*y(2)-(1/20)*y(2),y(1)^2-(1/20)*y(1),0];

y0=ones(3,1)*sqrt(3)/3;
tstar=1;
var=false;
ts=600;
while var==false
  
  k=tstar/ts;
  tt=linspace(0,1,ts+1);
  Tol=k^2;
				%Eulero implicito
  y=NaN(3,ts+1);
  y(:,1)=y0;
  t=0;
  for n=1:ts
    F=@(x) x-y(:,n)-k*f(t+k,x);
    JF=@(x) eye(3)-k*df(t+k,x);
    x0=y(:,n)+k*f(t,y(:,n));
    res=-JF(x0)\F(x0);
    while (norm(res,Inf)>Tol)
      x0+=res;
      res=-JF(x0)\F(x0);
    end
    x0+=res;
    y(:,n+1)=x0;
    t+=k;
  endfor
  
  mTilda1=sqrt(y(1,end)^2+y(2,end)^2+y(3,end)^2); 
  mTilda0=sqrt(y(1,1)^2+y(2,1)^2+y(3,1)^2); %mTilda0=1
  dist=abs(mTilda1-mTilda0)
  if (dist<2e-5)
    var=true;
    disp(ts)
  else
    var=false;
    ts+=20;
  endif

 % figure 
 % M=sqrt(y(1,:).^2+y(2,:).^2+y(3,:).^2);
 % plot(tt,M,'*');
 % hold on
endwhile


%% Verifica dell'implementazione (non richiesta)
counter=0;
tsrange=[10:10:100];
yRif=y;

for ts=tsrange
  k=tstar/ts;
  counter++;
  Tol=k^2;

  y=NaN(3,ts+1);
  y(:,1)=y0;
  t=0;
  for n=1:ts
    F=@(x) x-y(:,n)-k*f(t+k,x);
    JF=@(x) eye(3)-k*df(t+k,x);
    x0=y(:,n)+k*f(t,y(:,n));
    res=-JF(x0)\F(x0);
    while (norm(res,Inf)>Tol)
      x0+=res;
      res=-JF(x0)\F(x0);
    end
    x0+=res;
    y(:,n+1)=x0;
    t+=k;
  endfor

  err(counter)=norm(y(:,end)-yRif(:,end),Inf);

endfor
figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-1),'g')
toc
