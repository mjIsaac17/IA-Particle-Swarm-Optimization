# --------------------------------------|
# IA. Particle Swarm Optimization       |
# By: Isaac Montes Jiménez              |
# Created: Saturday October 12th, 2019  |
# Modified:                             |
# --------------------------------------|
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random as ran
import numpy as np
import math 

#Define the range of the values by the objective function 
#Ackley function
MAX_R_ACKLEY = 3
MIN_R_ACKLEY = -3
MAX_V = 0.8
ACKLEY_SOLUTION = 0

# C1 = the acceleration factor related to personal best
# C2 = the acceleration factor related to global best
C1 = 0.5
C2 = 1.5

# Store the best position found
V_GLOBAL_BEST = []

# Amount of particles
PARTICLES = 100

# Number of iterations
ITE = 300
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
    

#fn_llenarVectores()
#fn_plot(vX, vY, vZ, 'prueba')

class Particle:
    def __init__(self):
        self.vectorX = []
        self.vectorpBest = []
        self.vectorX.append(ran.uniform(MIN_R_ACKLEY, MAX_R_ACKLEY)) # current X position
        self.vectorX.append(ran.uniform(MIN_R_ACKLEY, MAX_R_ACKLEY)) # current Y position
        self.vectorpBest.append(self.vectorX[0]) # stores the position of the best solution found so far (by this particle)
        self.vectorpBest.append(self.vectorX[1])
        self.xFitness = fn_ackley_function(self.vectorX[0], self.vectorX[1])  # stores the current fitness of the particle
        self.pBestFitness = self.xFitness # stores the best fitness found by the particle
        self.velocity = ran.uniform((-1*MAX_V), MAX_V) # stores the gradient (direction) to move

    def setNextVelocity(self):
        r1 = ran.random()
        r2 = ran.random()
        self.velocity = self.velocity*0.4 + C1*r1*((self.vectorpBest[0] - self.vectorX[0]) + (self.vectorpBest[1] - self.vectorX[1])) + C2*r2*((V_GLOBAL_BEST[0] - self.vectorX[0]) + (V_GLOBAL_BEST[1] - self.vectorX[1]))

    def setNextPosition(self):
        self.vectorX[0] += self.velocity
        self.vectorX[1] += self.velocity

    def setFitness(self, newFitness):
        self.xFitness = newFitness
    
    def setBestFitness(self, newBestFitness):
        self.pBestFitness = newBestFitness


def main_PSO_ackley():
    particles = []
    #initialize the reference particle
    particles.append(Particle())
    bestAbsFitness = math.fabs(particles[0].xFitness) 
    V_GLOBAL_BEST.append(particles[0].vectorpBest[0])
    V_GLOBAL_BEST.append(particles[0].vectorpBest[1])
    #initialize the particles
    for i in range(1, PARTICLES):     
        particles.append(Particle())  
        print(i)
        if(math.fabs(particles[i].xFitness) < bestAbsFitness):
            bestAbsFitness = math.fabs(particles[i].xFitness)
            V_GLOBAL_BEST.append(particles[i].vectorpBest[0])
            V_GLOBAL_BEST.append(particles[i].vectorpBest[1])

    ite = 0
    while(ite < ITE):
        for i in range(PARTICLES):
            particles[i].setNextVelocity()
            particles[i].setNextPosition()
            particles[i].setFitness(fn_ackley_function(particles[i].vectorX[0], particles[i].vectorX[1]))
            # store the best fitness and position
            newAbsFitness = math.fabs(fn_ackley_function(particles[i].vectorX[0], particles[i].vectorX[1]))
            if(newAbsFitness < math.fabs(particles[i].pBestFitness)):
                particles[i].setBestFitness(fn_ackley_function(particles[i].vectorX[0], particles[i].vectorX[1]))
                particles[i].vectorpBest[0] = particles[i].vectorX[0]
                particles[i].vectorpBest[1] = particles[i].vectorX[1]
        ite += 1
    print("Best solution: ", fn_ackley_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1]))
    print("X: ", V_GLOBAL_BEST[0])
    print("Y: ", V_GLOBAL_BEST[1])

main_PSO_ackley()

#print("FN: ", fn_ackley_function(0.26663023740766967,1.3557306612252518))