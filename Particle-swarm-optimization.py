# --------------------------------------|
# IA. Particle Swarm Optimization       |
# By: Isaac Montes Jiménez              |
# Created: Saturday October 12th, 2019  |
# Modified:                             |
# --------------------------------------|
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import math 
"""
COUNT = 100
fig, ax = plt.subplots()
line, = ax.plot([], [], 'b.')
ax.set_ylim([-1.5, 1.5])
ax.set_xlim(0, COUNT)
xdata = []
ydata = [] 
def next():
   i = 0
   while i <= COUNT:
      i += 1
      yield i
def update(i):
   #xdata.append(i)
   xdata = i
   y = np.sin(i / 10.)
   #ydata.append(y)
   ydata = y
   line.set_data(xdata, ydata)
   return line,

#a = animation.FuncAnimation(fig, update, next, blit = False, interval = 10, repeat = False)
#plt.show()
"""
# objective function (ackley function)
def fn_ackley_function(x,y):
    f = -20*math.exp(-0.2*math.sqrt(0.5*(x**2+y**2))) - math.exp(0.5*(math.cos(2*math.pi*x) + math.cos(2*math.pi*y))) + math.e + 20
    return f

vX = []
vY = []
vZ = []
rango = 3
salto = 0.1

def fn_llenarVectores():
    x = -1 * rango
    while x < rango:
        y=-1 * rango
        while y < rango:
            vZ.append(fn_ackley_function(x,y))        
            vX.append(x)
            vY.append(y)
            y += salto
        x+= salto


#Plot the results
def fn_plot(array_x, array_y, array_energy, title):   
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.set_title(title)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z') 
    ax.plot_trisurf(array_x, array_y,array_energy, cmap=plt.cm.Spectral, antialiased=False)
    plt.show()
    

fn_llenarVectores()
fn_plot(vX, vY, vZ, 'prueba')