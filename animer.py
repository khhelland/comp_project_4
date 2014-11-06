from numpy import *
from matplotlib.pyplot import *
from sys import argv



filename = str(argv[1])

v = loadtxt(filename)
rand = zeros(len(v[:,0]))
v = column_stack((rand,v,rand))
x = linspace(0,1,11)


ion()
figure()
line, = plot(x,v[0,:])
xlim([0,1])
draw()

for w in v[:]:
    line.set_ydata(w)
    draw()
ioff()
show()
