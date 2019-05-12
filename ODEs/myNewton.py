import numpy as np

def Newton(f,df,x0,tol,maxiter):
    iter=0
    delta=np.linalg.solve(df(x0),-f(x0))
    while (np.linalg.norm(delta)>tol and iter<maxiter):
        iter+=1
        x0=x0+delta
        delta=np.linalg.solve(df(x0),-f(x0))
    x0=x0+delta
    
    if (iter>=maxiter):
        print('Massimo numero di iterazioni eseguito:', (iter))
        
    return x0

#Example
#def f(x):
#    return x**3-1
#def df(x):
#    return 3*x**2
#x0=-2
#tol=1e-5
#
#sol=Newton(f,df,x0,tol,200)