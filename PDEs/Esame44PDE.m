clear all; close all

m=150;
h=1/(m-1);
x=linspace(0,1,m)';

A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));

%% condizioni al bordo
A(1,1:2)=[-2,2]/h^2;
A(m,m-1:m)=[0,0];
%devo ancora imporre  che il termine di reazione sia con la prima componente nulla, in ULTIMA RIGA.

%termine reazione
g=@(u) [sin(u(1:m-1));0];
dg=@(u)[cos(u(1:m-1));0];

%%Ora ho y'=1/100A*y + sin(y);

%utilizzo punto medio implicito. Il problema è NON-lineare => Newton.
ts=1500;
tstar=1;
k=tstar/ts;
tol=k^2;
y0=5*(x.^2).*(x-1);

yn=y0;

F=@(y,yn) y - yn -k*(0.01*A*((y+yn)/2) + g((y+yn)/2));
J=@(y,yn) eye(length(y)) - k*0.5*0.01*A -k*0.5*spdiags(dg(y+yn),0,m,m);

for n=1:ts
	y=yn;
	res=-J(y,yn)\F(y,yn);
	
	while(norm(res,inf)>tol)
	y+=res;
	res=-J(y,yn)\F(y,yn);
	end
	y+=res;
	yn=y;	

	%plot(x,y)
	%axis([0,1,-1,1])	
	%pause(0.1)
end
yRef=y;

figure
plot(x,y)

counter=0;
tsrange=[100:100:500];

for ts=tsrange
	counter+=1;
	k=tstar/ts;
	tol=k^2;
	
	yn=y0;
	F=@(y,yn) y - yn -k*(0.01*A*((y+yn)/2) + g((y+yn)/2));
	J=@(y,yn) eye(length(y)) - k*0.5*0.01*A -k*0.5*spdiags(dg(y+yn),0,m,m);

	for n=1:ts
		y=yn;
		res=-J(y,yn)\F(y,yn);
		
		while(norm(res,inf)>tol)
		y+=res;
		res=-J(y,yn)\F(y,yn);
		end
		y+=res;
		yn=y;	
	end
		
	err(counter)=norm(y-yRef,inf);
end
figure
loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2),'r')
