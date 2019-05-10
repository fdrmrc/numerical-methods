clear all
close all

m = 100;
l = 1;
g = 9.81;
%accelerazione di gravità)
tstar = pi*sqrt(l/g);
k = tstar/m;
theta0 = pi/4;
y = [theta0; 0];
% inizializzo le funzioni per pendolo non linearizzato
f = @(y) [y(2); -(g/l)*sin(y(1))];
jf = @(y) [0 1 ; 0 -(g/l)*cos(y(1))];
% METODO DEI TRAPEZI
F = @(y,yn) y-yn-k/2*(f(yn)+f(y));
JF = @(y) eye(length(y)) - k/2*jf(y);


% inizializzo le funzioni per il pendolo linearizzato
A = [0 1 ; -g/l, 0];
y2 = y(:,1);
I = eye(length(y2));

tol = k^2/100;
t = 0;

for i = 1:m
    % risoluzione non linearizzato: applico Newton
    t =  t+k;
    % guess iniziale
    y(:,i+1) = y(:,i);
    
    delta = -JF(y(:,i+1))\F(y(:,i+1),y(:,i));
    while norm(delta)>tol
        y(:,i+1) = y(:,i+1)+delta;
        delta = -JF(y(:,i+1))\F(y(:,i+1),y(:,i));
    end 
    % risoluzione linearizzato
    y2(:,i+1) = (I-k/2*A)\((I+k/2*A)*y2(:,i));
end       
       
% Trovare il numero minimo di passi affichè la soluzione attraverso Eulero
% disti da theta(tstar) meno di 10^-1

% Guess iniziali
tol = 1;
meul = 10;

while tol > 1e-2
    keul = tstar/meul;
    yeul = [theta0;0];
    for j = 1:meul
        % Applico Eulero esplicito
        yeul(:,j+1) = yeul(:,j)+keul*A*yeul(:,j);
    end
    % Trovo nuova tolleranza e aggiorno meul
    tol = abs(yeul(1,end)-y2(1,end));
    meul = meul+1;
 end
meul = meul-1;
   
subplot(1,2,1)
plot(0:k:tstar,y(1,:),'x',0:k:tstar,y2(1,:),'o')
subplot(1,2,2)
plot(0:k:tstar,y2(1,:),'o',0:keul:tstar,yeul(1,:),'.')
