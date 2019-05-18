import numpy as np
import matplotlib.pyplot as plt
import math as m

#Trapz rule is NOT L-stable! (While Backward Euler is L-Stable)
#Test equation: y'(t)=-2000*(y-cos(t)), y(0)=0

def Fun(t,x):
    return -2000*(x-m.cos(t))

tf=1.5
t0=0
ts=50
k=(tf-t0)/ts
t=np.linspace(t0,tf,ts+1)
y=np.zeros((ts+1)) #trapz
yE=np.zeros((ts+1)) #Backward Euler
y0=0
y[0]=y0
yE[0]=0

for i in range (0,ts):
    y[i+1]=(y[i]+1000*k*(-y[i]+m.cos(t[i])+m.cos(t[i+1])))/(1+k*1000)
    yE[i+1]=(y[i]+2000*k*m.cos(t[i+1]))/(1+k*2000)

plt.plot(t,y,marker='o',label='Trapezi')        
plt.plot(t,yE,marker='d',label='Backward Euler')
plt.legend()
plt.show()