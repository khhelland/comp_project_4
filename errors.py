import numpy as np
import matplotlib.pyplot as plt
import os,sys
from math import pi

def sumx(n):
    a = 2*((np.sin(n*pi) - n*pi)/(n*pi)**2)*np.sin(n*pi*x)
    return a
def sumt(n,t):
    return np.exp(-t*(n*pi)**2)


error = os.system('make')
if error:
    print 'make funket ikke'
    sys.exit(1)

run = './main.x'
error = os.system(run)
if error:
    print 'runtime error'
    sys.exit(1)

x = np.linspace(0,1,11)[1:-1]
t = np.array((0.02,0.5))

v = np.zeros((9,2))


for m in range(1,300):
    v[:,0] += sumx(m)*sumt(m,t[0])
    v[:,1] += sumx(m)*sumt(m,t[1])

FE = np.loadtxt('ForwardEuler.dat').T
BE = np.loadtxt('BackwardEuler.dat').T
CN = np.loadtxt('CrankNicholson.dat').T

dFE = np.abs((FE)/v-1)
dBE = np.abs((BE)/v-1)
dCN = np.abs((CN)/v-1)

mdfe0 = max(dFE[0])
mdfe1 = max(dFE[1])

mdbe0 = max(dBE[0])
mdbe1 = max(dBE[1])

mdcn0 = max(dCN[0])
mdcn1 = max(dCN[1])


print 'maximal relative errors for FE,BE and CN at t = 0.02:'
print mdfe0, '\t', mdbe0, '\t',mdcn0 

print 'maximal relative errors for FE,BE and CN at t = 0.5:'
print mdfe1, '\t', mdbe1, '\t',mdcn1 


plt.figure()
plt.plot(x,v[:,0],'-x',label="sol")
plt.hold(1)
plt.plot(x,FE[:,0],'-x', label="FE")
plt.plot(x,BE[:,0],'-x',label="BE")
plt.plot(x,CN[:,0],'-x',label="CN")
plt.legend(loc=4)

plt.figure()
plt.plot(x,v[:,1],'-x',label="sol")
plt.hold(1)
plt.plot(x,FE[:,1],'-x',label="FE")
plt.plot(x,BE[:,1],'-x',label="BE")
plt.plot(x,CN[:,1],'-x',label="CN")
plt.legend(loc=4)

plt.show()
