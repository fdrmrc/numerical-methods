import numpy as np
import matplotlib.pyplot as plt
import math

A=np.array([[1,0],[-1,0]])

def f(t,y):    
    return np.matmul(A,y)+[0,2*math.cos(t)]


ts=10000
tf=1
k=tf/ts
y0=np.array([1,0])
y=np.zeros((2,ts+1))
y[:,0]=y0
t=np.linspace(0,tf,ts+1)

for i in range (0,ts):
    y[:,i+1]=np.linalg.solve(np.eye(2)-k*A,y[:,i]+k*np.array([0,2*math.cos(t[i+1])]))

plt.plot(t,y[0,:],linestyle=':',label='sol numerica')
plt.plot(t,np.exp(t),label='analitica')
plt.axis([0,tf,np.min(y[0,:]),np.max(y[0,:])])
plt.legend()
plt.show()

Err=np.linalg.norm(y[0,:]-np.exp(t))
print('Errore wrt soluzione analitica=\n',Err)