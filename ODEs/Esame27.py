import numpy as np
import matplotlib.pyplot as plt
import math
from myNewton import Newton
def fun(x):
    return np.array([-2*x[0]*x[1],x[0]**2+x[2]**2-x[1]**2 -1,-2*(x[1]+x[0])*x[2]])

def dfun(x):
    return np.array([[-2*x[1],-2*x[0],0],[2*x[0],-2*x[1],2*x[2]],[-2*x[2],-2*x[2],-2*(x[0]+x[1])]])
    
tf=1
tsrif=2000
ts=tsrif
k=tf/ts
y0=np.array([1,1,1])
y=np.zeros((len(y0),ts+1))
y[:,0]=y0
tol=0.01*k**2 
for i in range (0,ts):
    yi=y[:,i]
    def F(x):
        return np.array(x-yi-0.5*k*(fun(yi)+fun(x)))
    def JF(x):
        return np.array(np.eye(3)-0.5*k*dfun(x))
    x0=yi
    y[:,i+1]=Newton(F,JF,x0,tol,100)

t=np.linspace(0,tf,ts+1)
#plt.plot(t,y[0,:])
#plt.show()
yRif=y

err=[]
tsrange=np.arange(10,510,10)
for ts in tsrange:
    k=tf/ts
    y=np.zeros((len(y0),ts+1))
    y[:,0]=y0
    tol=0.01*k**2 
    for i in range (0,ts):
        yi=y[:,i]
        def F(x):
            return np.array(x-yi-0.5*k*(fun(yi)+fun(x)))
        def JF(x):
            return np.array(np.eye(3)-0.5*k*dfun(x))
        x0=yi
        y[:,i+1]=Newton(F,JF,x0,tol,100)
    err.append(np.linalg.norm(yRif[:,tsrif]-y[:,ts]))

plt.loglog(tsrange,err,'*',label='ordine effettivo')
plt.loglog(tsrange,err[-1]*np.power((tsrange/tsrange[-1]),-2),label='$\mathcal{O}(h^2)$')
plt.legend()
plt.show()
    