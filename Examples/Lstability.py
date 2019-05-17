import numpy as np
import matplotlib.pyplot as plt
import math as m

def Fun(t,x):
    return -2000*(x-m.cos(t))

tf=3/2
t0=0
ts=50
k=(tf-t0)/ts
y=np.zeros((ts+1))
y0=0
y[0]=y0
t=np.linspace(t0,tf,ts+1)


for i in range (0,ts):
    y[i+1]=(y[i]+1000*k*(y[i]+m.cos(t[i])+m.cos(t[i+1])))/(1-k*1000)

plt.plot(t,y,marker='o')
plt.show()