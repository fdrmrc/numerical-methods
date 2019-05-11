import math
import numpy as np
import matplotlib.pyplot as plt

def fun(x):
    return np.array([-2*x[0]*x[1],x[0]**2+x[2]**2-x[1]**2 -1,-2*(x[1]+x[0])*x[2]])

def dfun(x):
    return np.array([[-2*x[1],-2*x[0],0],[2*x[0],-2*x[1],2*x[2]],[-2*x[2],-2*x[2],-2*(x[0]+x[1])]])

y0=np.array([1,2,15])
tsrif=1000
ts=tsrif
tf=1
k=tf/ts
y=np.zeros((3,ts+1))
y[:,0]=y0
Tol=0.01*k**3
for i in range (0,ts):
    yi=y[:,i]
    def F(x):
        return np.array(x-yi-0.5*k*(fun(yi)+fun(x)))
    
    def JF(x):
        return np.array(np.eye(3)-0.5*k*dfun(x))
    
    x0=yi
    delta=np.linalg.solve(JF(x0),-F(x0))
    while (np.linalg.norm(delta)>Tol):
        x0=x0+delta
        delta=np.linalg.solve(JF(x0),-F(x0))
    x0=x0+delta
    y[:,i+1]=x0

#t=np.linspace(0,tf,ts+1)
#plt.plot(t,y[0,:])
#plt.show()
    
yRif=y #soluzione di riferimento

#Mostro ordine di convergenza
tsrange=np.arange(50,350,50)
err=[]
for ts in tsrange:
    k=tf/ts
    y=np.zeros((3,ts+1))
    y[:,0]=y0
    Tol=0.01*k**3
    for i in range (0,ts):
        yi=y[:,i]
        def F(x):
            return np.array(x-yi-0.5*k*(fun(yi)+fun(x)))
        def JF(x):
            return np.array(np.eye(3)-0.5*k*dfun(x))
        x0=yi
        delta=np.linalg.solve(JF(x0),-F(x0))
        while (np.linalg.norm(delta)>Tol):
            x0=x0+delta
            delta=np.linalg.solve(JF(x0),-F(x0))
        
        x0=x0+delta
        y[:,i+1]=x0
    err.append(np.linalg.norm(y[:,ts]-yRif[:,tsrif]))

plt.loglog(tsrange,err,'*',label='ordine di convergenza effettivo')
plt.loglog(tsrange,err[-1]*1/np.power((tsrange/tsrange[-1]),2),marker='>',markersize=2,label='ordine di convergenza atteso=$\mathcal{O}(h^2)$')
plt.legend()
plt.xlabel('timesteps')
plt.ylabel('errore')
plt.show()
