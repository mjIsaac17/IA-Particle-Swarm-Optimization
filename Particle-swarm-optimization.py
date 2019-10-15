# IA - Particle Swarm Optimization 
# By: Isaac Montes Jiménez
# Created: 10/13/2019
# Modified: 10/15/2019

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D 
import numpy as np
import random as ran
import math 

#Define the range of the values by the objective function 
#Ackley function
MAX_R_ACKLEY = 4
MIN_R_ACKLEY = -4

#Beale function
MAX_R_BEALE = 4.5
MIN_R_BEALE = -4.5

# Max initial particle velocity
MAX_V = 2

# C1 = the acceleration factor related to personal best
# C2 = the acceleration factor related to global best
C1 = 1
C2 = 2

# Store the best position found
V_GLOBAL_BEST = []

# Amount of particles
PARTICLES = 200

# Number of iterations
ITE = 50

# Initial weight
W = [0.9]
W_REDUCTION = 0.01 # Reduction factor to each iteration 

# Flag to know if the model of the objective function is plotted
# 0 = OFF
# 1 = ON
MODEL = 0

# Value between points of the model, lower = more quality
MODEL_QUALITY = 0.3

# Speed of the iterations, lower = faster
SPEED = 100

# Flag to know the objective function
# 1 = ACKLEY
# 2 = BEALE
OBJ_FUNCTION = 2

# Objective function (Ackley)
def fn_ackley_function(x,y):
    f = -20*math.exp(-0.2*math.sqrt(0.5*(x**2+y**2))) - math.exp(0.5*(math.cos(2*math.pi*x) + math.cos(2*math.pi*y))) + math.e + 20
    return f

# Objective function (Beale)
def fn_beale_function(x,y):
    f = (1.5 - x + x*y)**2 + (2.25 - x + (x*y**2))**2 + (2.625 - x + (x*y**3))**2
    return f

class Particle:
    def __init__(self):
        self.vectorX = []
        self.vectorpBest = []
        self.vectorX.append(ran.uniform(MIN_R_ACKLEY, MAX_R_ACKLEY)) # current X position
        self.vectorX.append(ran.uniform(MIN_R_ACKLEY, MAX_R_ACKLEY)) # current Y position
        self.vectorpBest.append(self.vectorX[0]) # stores the position of the best solution found so far (by this particle)
        self.vectorpBest.append(self.vectorX[1])
        if(OBJ_FUNCTION == 1):
            self.xFitness = fn_ackley_function(self.vectorX[0], self.vectorX[1])  # stores the current fitness of the particle
        elif(OBJ_FUNCTION == 2):
            self.xFitness = fn_beale_function(self.vectorX[0], self.vectorX[1])
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


def fn_createModelAckley():
    vX = []
    vY = []
    vZ = []
    x = -1 * MAX_R_ACKLEY
    while x < MAX_R_ACKLEY:
        y=-1 * MAX_R_ACKLEY
        while y < MAX_R_ACKLEY:
            vZ.append(fn_ackley_function(x,y))        
            vX.append(x)
            vY.append(y)
            y += MODEL_QUALITY
        x+= MODEL_QUALITY
    return vX, vY, vZ

def fn_createModelBeale():
    vX = []
    vY = []
    vZ = []
    x = -1 * MAX_R_BEALE
    while x < MAX_R_BEALE:
        y=-1 * MAX_R_BEALE
        while y < MAX_R_BEALE:
            vZ.append(fn_beale_function(x,y))        
            vX.append(x)
            vY.append(y)
            y += MODEL_QUALITY
        x+= MODEL_QUALITY
    return vX, vY, vZ


def main_PSO_3D():
    particles = []
    modelVX = []
    modelVY = []
    modelVZ = []
    modelInfo = []
    particles.append(Particle())
    bestAbsFitness = math.fabs(particles[0].xFitness) 
    V_GLOBAL_BEST.append(particles[0].vectorpBest[0])
    V_GLOBAL_BEST.append(particles[0].vectorpBest[1])
    MAX_PLOT_R = 0
    MIN_PLOT_R = 0

    if(OBJ_FUNCTION == 1): #Ackley function
        modelInfo = ['Ackley function', 'Solución: 0', 'x = 0', 'y = 0']
        V_GLOBAL_BEST.append(fn_ackley_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1]))
        MAX_PLOT_R = MAX_R_ACKLEY
        MIN_PLOT_R = MIN_R_ACKLEY
        #initialize the particles
        for i in range(1, PARTICLES):     
            particles.append(Particle())  
            if(math.fabs(particles[i].xFitness) < bestAbsFitness):
                bestAbsFitness = math.fabs(particles[i].xFitness)
                V_GLOBAL_BEST[0] = particles[i].vectorpBest[0]
                V_GLOBAL_BEST[1] = particles[i].vectorpBest[1]
                V_GLOBAL_BEST[2] = fn_ackley_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1]) 

    elif(OBJ_FUNCTION == 2): #Beale function
        modelInfo = ['Beale function', 'Solución: 0', 'x = 3', 'y = 0.5']
        V_GLOBAL_BEST.append(fn_beale_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1]))
        MAX_PLOT_R = MAX_R_BEALE
        MIN_PLOT_R = MIN_R_BEALE
        #initialize the particles
        for i in range(1, PARTICLES):     
            particles.append(Particle())  
            if(math.fabs(particles[i].xFitness) < bestAbsFitness):
                bestAbsFitness = math.fabs(particles[i].xFitness)
                V_GLOBAL_BEST[0] = particles[i].vectorpBest[0]
                V_GLOBAL_BEST[1] = particles[i].vectorpBest[1]
                V_GLOBAL_BEST[2] = fn_beale_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1]) 

    # Fill the first position of the particles in the animation
    x = []
    y = []
    z = []
    for i in range(PARTICLES):
        x.append(particles[i].vectorX[0])
        y.append(particles[i].vectorX[1])
        z.append(particles[i].xFitness)

    # Create the figure to animate
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    sct, = ax.plot([], [], [], "ok", markersize=4)
    ax.grid(True, linestyle = '-', color = '0.75')
    scale = 1
    ax.set_xlim([scale * MIN_PLOT_R, scale * MAX_PLOT_R])
    ax.set_ylim([scale * MIN_PLOT_R, scale * MAX_PLOT_R])
    ax.set_zlim([-1, 175000])    
    if(MODEL):
        if(OBJ_FUNCTION == 1):
            modelVX, modelVY, modelVZ = fn_createModelAckley()
        elif (OBJ_FUNCTION == 2):
            modelVX, modelVY, modelVZ = fn_createModelBeale()
        ax.plot_trisurf(modelVX, modelVY, modelVZ, cmap=plt.cm.Spectral, antialiased=False)
    
    def _update_plot(j): 
        plt.title("Modelo: %s / %s / %s , %s \n\nIteraciones: %d / Mejor Fitness: %f / x = %f , y = %f" 
            %(modelInfo[0], modelInfo[1], modelInfo[2], modelInfo[3], j, V_GLOBAL_BEST[2], V_GLOBAL_BEST[0], V_GLOBAL_BEST[1])) 
        if(j == ITE):
            anim._stop()
        partX = []
        partY = []
        partZ = []
        absBestFitness = math.fabs(V_GLOBAL_BEST[2])

        if(OBJ_FUNCTION == 1):
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
                
                # Check the best particle to update the reference of the best one
                if(newAbsFitness < absBestFitness):
                    V_GLOBAL_BEST[0] = particles[i].vectorpBest[0]
                    V_GLOBAL_BEST[1] = particles[i].vectorpBest[1]
                    V_GLOBAL_BEST[2] = fn_ackley_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1]) 
                    absBestFitness = math.fabs(V_GLOBAL_BEST[2])
                partX.append(particles[i].vectorX[0])
                partY.append(particles[i].vectorX[1])
                partZ.append(particles[i].xFitness)

        elif (OBJ_FUNCTION == 2):
            for i in range(PARTICLES):
                particles[i].setNextVelocity()
                particles[i].setNextPosition()
                particles[i].setFitness(fn_beale_function(particles[i].vectorX[0], particles[i].vectorX[1]))
                # store the best fitness and position
                newAbsFitness = math.fabs(fn_beale_function(particles[i].vectorX[0], particles[i].vectorX[1]))
                if(newAbsFitness < math.fabs(particles[i].pBestFitness)):
                    particles[i].setBestFitness(fn_beale_function(particles[i].vectorX[0], particles[i].vectorX[1]))
                    particles[i].vectorpBest[0] = particles[i].vectorX[0]
                    particles[i].vectorpBest[1] = particles[i].vectorX[1]
                
                # Check the best particle to update the reference of the best one
                if(newAbsFitness < absBestFitness):
                    V_GLOBAL_BEST[0] = particles[i].vectorpBest[0]
                    V_GLOBAL_BEST[1] = particles[i].vectorpBest[1]
                    V_GLOBAL_BEST[2] = fn_beale_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1]) 
                    absBestFitness = math.fabs(V_GLOBAL_BEST[2])
                partX.append(particles[i].vectorX[0])
                partY.append(particles[i].vectorX[1])
                partZ.append(particles[i].xFitness)
        sct.set_data(partX, partY)
        sct.set_3d_properties(partZ)
        W[0] -= W_REDUCTION

    sct.set_data(x, y)
    sct.set_3d_properties(z)
    anim = animation.FuncAnimation(fig, _update_plot, interval = SPEED)
    plt.show()


def main_PSO_2D():
    particles = []
    particles.append(Particle())
    bestAbsFitness = math.fabs(particles[0].xFitness) 
    V_GLOBAL_BEST.append(particles[0].vectorpBest[0])
    V_GLOBAL_BEST.append(particles[0].vectorpBest[1])
    modelInfo = []
    if(OBJ_FUNCTION == 1):
        modelInfo = ['Ackley function', 'Solución: 0', 'x = 0', 'y = 0']
        V_GLOBAL_BEST.append(fn_ackley_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1]))
        MAX_PLOT_R = MAX_R_ACKLEY
        MIN_PLOT_R = MIN_R_ACKLEY
        #initialize the particles
        for i in range(1, PARTICLES):     
            particles.append(Particle())  
            if(math.fabs(particles[i].xFitness) < bestAbsFitness):
                bestAbsFitness = math.fabs(particles[i].xFitness)
                V_GLOBAL_BEST[0] = particles[i].vectorpBest[0]
                V_GLOBAL_BEST[1] = particles[i].vectorpBest[1]
                V_GLOBAL_BEST[2] = fn_ackley_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1])
    elif(OBJ_FUNCTION == 2):
        modelInfo = ['Beale function', 'Solución: 0', 'x = 3', 'y = 0.5']
        V_GLOBAL_BEST.append(fn_beale_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1]))
        MAX_PLOT_R = MAX_R_BEALE
        MIN_PLOT_R = MIN_R_BEALE
        #initialize the particles
        for i in range(1, PARTICLES):     
            particles.append(Particle())  
            if(math.fabs(particles[i].xFitness) < bestAbsFitness):
                bestAbsFitness = math.fabs(particles[i].xFitness)
                V_GLOBAL_BEST[0] = particles[i].vectorpBest[0]
                V_GLOBAL_BEST[1] = particles[i].vectorpBest[1]
                V_GLOBAL_BEST[2] = fn_beale_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1])
    x = []
    y = []
    for i in range(PARTICLES):
        x.append(particles[i].vectorX[0])
        y.append(particles[i].vectorX[1])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid(True, linestyle = '-', color = '0.75')
    scale = 2
    ax.set_xlim([scale * MIN_PLOT_R, scale * MAX_PLOT_R])
    ax.set_ylim([scale * MIN_PLOT_R, scale * MAX_PLOT_R])
    sct, = ax.plot([], [], "or", markersize=3) # o = draw dots, r = red, red dots

    def update_plot(j): 
        partX = []
        partY = []
        plt.title("Modelo: %s / %s / %s , %s \n\nIteraciones: %d / Mejor Fitness: %f / x = %f , y = %f" 
            %(modelInfo[0], modelInfo[1], modelInfo[2], modelInfo[3], j, V_GLOBAL_BEST[2], V_GLOBAL_BEST[0], V_GLOBAL_BEST[1])) 
        if(j == ITE):
            anim._stop()
        absBestFitness = math.fabs(V_GLOBAL_BEST[2])
        if(OBJ_FUNCTION == 1):
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

                # Check the best particle to update the reference of the best one
                if(newAbsFitness < absBestFitness):
                    V_GLOBAL_BEST[0] = particles[i].vectorpBest[0]
                    V_GLOBAL_BEST[1] = particles[i].vectorpBest[1]
                    V_GLOBAL_BEST[2] = fn_ackley_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1]) 
                    absBestFitness = math.fabs(V_GLOBAL_BEST[2])
                partX.append(particles[i].vectorX[0])
                partY.append(particles[i].vectorX[1])
        elif(OBJ_FUNCTION == 2):
            for i in range(PARTICLES):
                particles[i].setNextVelocity()
                particles[i].setNextPosition()
                particles[i].setFitness(fn_beale_function(particles[i].vectorX[0], particles[i].vectorX[1]))
                # store the best fitness and position
                newAbsFitness = math.fabs(fn_beale_function(particles[i].vectorX[0], particles[i].vectorX[1]))
                if(newAbsFitness < math.fabs(particles[i].pBestFitness)):
                    particles[i].setBestFitness(fn_beale_function(particles[i].vectorX[0], particles[i].vectorX[1]))
                    particles[i].vectorpBest[0] = particles[i].vectorX[0]
                    particles[i].vectorpBest[1] = particles[i].vectorX[1]

                # Check the best particle to update the reference of the best one
                if(newAbsFitness < absBestFitness):
                    V_GLOBAL_BEST[0] = particles[i].vectorpBest[0]
                    V_GLOBAL_BEST[1] = particles[i].vectorpBest[1]
                    V_GLOBAL_BEST[2] = fn_beale_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1]) 
                    absBestFitness = math.fabs(V_GLOBAL_BEST[2])
                partX.append(particles[i].vectorX[0])
                partY.append(particles[i].vectorX[1])
        sct.set_data(partX, partY)
        W[0] -= W_REDUCTION

    sct.set_data(x, y)
    anim = animation.FuncAnimation(fig, update_plot,interval = SPEED)
    plt.show()

main_PSO_2D()
#main_PSO_3D()

