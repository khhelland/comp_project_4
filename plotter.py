import numpy as np
import matplotlib.pyplot as plt
import os,sys
from math import pi


def sumx(n):
    a = ((np.sin(n*pi) - n*pi)/(n*pi)**2)*np.sin(n*pi*x)
    return a
def sumt(n,t):
    return np.exp(-t*(n*pi)**2)

x = np.linspace(0,1,100)
t = np.array((0.02,0.5))

v = np.zeros((100,2))


for m in range(1,300):
    v[:,0] += sumx(m)*sumt(m,t[0])
    v[:,1] += sumx(m)*sumt(m,t[1])

plt.figure()
plt.plot(x,v[:,0], label = "time = %g"%t[0])
plt.hold(1)
plt.plot(x,v[:,1], label = "time = %g"%t[1])

plt.xlabel('x')
plt.ylabel('v')
plt.legend(loc = 4)
plt.show()
