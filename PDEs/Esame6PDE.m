clear all
close all
m = 500;
x = linspace(0,1,m)';
h = 1/(m-1);
A = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
A(m,m-1:m)=[2,-2]/h^2;
A(1,1:2)=[0,0];
b=@(t) [0;-sin(t)*[-x(2:m).^2 + 2*x(2:m)] + 2*cos(t)];

y0 = -x.^2 + 2*x;
tsrange = 10:10:100;
counter = 0;
for ts = tsrange
    counter = counter+1;
    k = 1/ts;
    y = y0;
    t = 0;
    P = phi1m(k*A);
    for n = 1:ts
      y = y+k*P*(A*y+b(t)); % eulero esponenziale
      t = t+k;
    end
    err(counter) = norm(y-cos(t)*(-x.^2+2*x),inf);
end

loglog(tsrange,err,'*',...
    tsrange,err(end)*(tsrange/tsrange(end)).^(-1))
title('Error plot')
xlabel('timesteps')
ylabel('Inf norm of the error')
