import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D 
import numpy as np
from numpy.random import normal as normal
import random as ran
import math 

#Define the range of the values by the objective function 
#Ackley function
MAX_R_ACKLEY = 3
MIN_R_ACKLEY = -3
MAX_V = 2
ACKLEY_SOLUTION = 0

# C1 = the acceleration factor related to personal best
# C2 = the acceleration factor related to global best
C1 = 1
C2 = 2

# Store the best position found
V_GLOBAL_BEST = []

# Amount of particles
PARTICLES = 100

# Number of iterations
ITE = 50
W = [0.9]
particles = []



def fn_ackley_function(x,y):
    f = -20*math.exp(-0.2*math.sqrt(0.5*(x**2+y**2))) - math.exp(0.5*(math.cos(2*math.pi*x) + math.cos(2*math.pi*y))) + math.e + 20
    return f

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
        self.velocityX = ran.uniform((-1*MAX_V), MAX_V) # stores the gradient (direction) to move
        self.velocityY = ran.uniform((-1*MAX_V), MAX_V) 

    def setNextVelocity(self):
        r1 = ran.random()
        r2 = ran.random()
        cognitivoX = C1*r1*(self.vectorpBest[0] - self.vectorX[0])
        cognitivoY = C1*r1*(self.vectorpBest[1] - self.vectorX[1])
        socialX = C2*r2*(V_GLOBAL_BEST[0] - self.vectorX[0])
        socialY = C2*r2*(V_GLOBAL_BEST[1] - self.vectorX[1])
        self.velocityX = W[0] * self.velocityX + (cognitivoX + socialX)
        self.velocityY = W[0] * self.velocityY + (cognitivoY + socialY)

    def setNextPosition(self):
        self.vectorX[0] += self.velocityX
        self.vectorX[1] += self.velocityY

    def setFitness(self, newFitness):
        self.xFitness = newFitness
    
    def setBestFitness(self, newBestFitness):
        self.pBestFitness = newBestFitness


vX = []
vY = []
vZ = []
rango = 3
salto = 0.2
MODELO = 1
band = 0

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
    band = 1

#Plot the results
def fn_plot(array_x, array_y, array_energy, title):   
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.set_title(title)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z') 
    ax.plot_trisurf(array_x, array_y,array_energy, cmap=plt.cm.Spectral, antialiased=False)


### INICIO
def ackley():
    particles.append(Particle())
    bestAbsFitness = math.fabs(particles[0].xFitness) 
    V_GLOBAL_BEST.append(particles[0].vectorpBest[0])
    V_GLOBAL_BEST.append(particles[0].vectorpBest[1])
    V_GLOBAL_BEST.append(fn_ackley_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1]))

    #initialize the particles
    for i in range(1, PARTICLES):     
        particles.append(Particle())  
        if(math.fabs(particles[i].xFitness) < bestAbsFitness):
            bestAbsFitness = math.fabs(particles[i].xFitness)
            V_GLOBAL_BEST[0] = particles[i].vectorpBest[0]
            V_GLOBAL_BEST[1] = particles[i].vectorpBest[1]
            V_GLOBAL_BEST[2] = fn_ackley_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1]) 
    x = []
    y = []
    z = []
    for i in range(PARTICLES):
        x.append(particles[i].vectorX[0])
        y.append(particles[i].vectorX[1])
        z.append(particles[i].xFitness)
    #ax = fig.add_subplot(111)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    sct, = ax.plot([], [], [], "o", markersize=2)
    ax.grid(True, linestyle = '-', color = '0.75')
    ax.set_xlim([-5, 5])
    ax.set_ylim([-5, 5])
    ax.set_zlim([-1, 8])

        
    if(MODELO):
        fn_llenarVectores()
        ax.plot_trisurf(vX, vY,vZ, cmap=plt.cm.Spectral, antialiased=False)
    #fn_plot(vX, vY, vZ, "i")
    #scat = plt.scatter(x, y, z)
    #scat.set_alpha(0.8)
    def _update_plot(j): 
        plt.title("Iteraciones: %d / Mejor Fitness: %f" %(j, V_GLOBAL_BEST[2])) 
        if(j == ITE):
            anim._stop()
        partX = []
        partY = []
        partZ = []
        absBestFitness = math.fabs(V_GLOBAL_BEST[2])
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
            if(newAbsFitness < absBestFitness):
                V_GLOBAL_BEST[0] = particles[i].vectorpBest[0]
                V_GLOBAL_BEST[1] = particles[i].vectorpBest[1]
                V_GLOBAL_BEST[2] = fn_ackley_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1]) 
                absBestFitness = math.fabs(V_GLOBAL_BEST[2])
            partX.append(particles[i].vectorX[0])
            partY.append(particles[i].vectorX[1])
            partZ.append(particles[i].xFitness)
        sct.set_data(partX, partY)
        sct.set_3d_properties(partZ)
        
        W[0] -= 0.01
    #scat.set_offsets(part)
    anim = animation.FuncAnimation(fig, _update_plot, interval = 80)

    plt.show()

ackley()
def grafica():

    nfr = 50 # Number of frames
    fps = 10 # Frame per sec
    xs = []
    ys = []
    zs = []
    ss = np.arange(1,nfr,0.5)
    for s in ss:
        xs.append(normal(0,s,200))
        ys.append(normal(50,s,200))
        zs.append(normal(50,s,200))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    sct, = ax.plot([], [], [], "o", markersize=2)
    def update(ifrm, xa, ya, za):
        sct.set_data(xa[ifrm], ya[ifrm])
        sct.set_3d_properties(za[ifrm])
    print(nfr)
    ax.set_xlim(0,100)
    ax.set_ylim(0,100)
    ax.set_zlim(0,100)
    ani = animation.FuncAnimation(fig, update, nfr, fargs=(xs,ys,zs), interval=1000/fps)
    plt.show()

