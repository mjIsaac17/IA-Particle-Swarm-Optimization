import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
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
PARTICLES = 10

# Number of iterations
ITE = 200
W = [0.9]
particles = []
def _update_plot(j, fig, scat): ##FALTA ACTUALIZAR LA MEJOR POSICION
    ite = 0
    part = []
    #while(ite < ITE):
    absBestFitness = math.fabs(V_GLOBAL_BEST[2])
    for i in range(PARTICLES):
        particles[i].setNextVelocity()
        particles[i].setNextPosition()
        #print("particle %d, posX: %f, posY: %f, velX: %f, velY: %f" %(i, particles[i].vectorX[0], particles[i].vectorX[1], particles[i].velocityX, particles[i].velocityY))
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
        part.append([particles[i].vectorX[0], particles[i].vectorX[1]])
        #ite += 1
    #print("BEST X: %f, BEST Y: %f" %{V_GLOBAL_BEST[0], V_GLOBAL_BEST[1]})
    #print("______________________________________")
    W[0] -= 0.01
    scat.set_offsets(part)
    #print('Frames: %d', %i) # For < Python3 
    #print('Frames: {}'.format(j)) # For Python3
    
    return scat


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
        #print("BEST X: ", V_GLOBAL_BEST[0], "BEST Y: ", V_GLOBAL_BEST[1])
        r1 = ran.random()
        r2 = ran.random()
        cognitivoX = C1*r1*(self.vectorpBest[0] - self.vectorX[0])
        cognitivoY = C1*r1*(self.vectorpBest[1] - self.vectorX[1])
        socialX = C2*r2*(V_GLOBAL_BEST[0] - self.vectorX[0])
        socialY = C2*r2*(V_GLOBAL_BEST[1] - self.vectorX[1])
        #print("COG: R1(%f)*(PBX(%f) + PBY(%f) - PX(%f) - PX(%f)) = %f" %(r1, self.vectorpBest[0], self.vectorpBest[1] , self.vectorX[0], self.vectorX[1], cognitivo))
        #print("SOC: R2(%f)*(GBX(%f) + GBY(%f) - PX(%f) - PX(%f)) = %f" %(r2, V_GLOBAL_BEST[0], V_GLOBAL_BEST[1] , self.vectorX[0], self.vectorX[1], social))
        #print("NEW VEL: VEL(%f)*0.5 + COG(%f) + SOC(%f) = %f / NEWX: %f, NEWY:%f" %(self.velocity, cognitivo, social, self.velocity*0.5+cognitivo+social, self.vectorX[0]+self.velocity*0.5+cognitivo+social,self.vectorX[1]+self.velocity*0.5+cognitivo+social))
        #self.velocity = self.velocity*0.4 + C1*r1*((self.vectorpBest[0] - self.vectorX[0]) - (self.vectorpBest[1] - self.vectorX[1])) + C2*r2*((V_GLOBAL_BEST[0] - self.vectorX[0]) - (V_GLOBAL_BEST[1] - self.vectorX[1]))
        self.velocityX = W[0] * self.velocityX + (cognitivoX + socialX)
        self.velocityY = W[0] * self.velocityY + (cognitivoY + socialY)
    def setNextPosition(self):
        self.vectorX[0] += self.velocityX
        self.vectorX[1] += self.velocityY

    def setFitness(self, newFitness):
        self.xFitness = newFitness
    
    def setBestFitness(self, newBestFitness):
        self.pBestFitness = newBestFitness


def main_PSO_ackley():
    
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
            V_GLOBAL_BEST[0] = particles[i].vectorpBest[0]
            V_GLOBAL_BEST[1] = particles[i].vectorpBest[1]
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



particles.append(Particle())
bestAbsFitness = math.fabs(particles[0].xFitness) 
V_GLOBAL_BEST.append(particles[0].vectorpBest[0])
V_GLOBAL_BEST.append(particles[0].vectorpBest[1])
V_GLOBAL_BEST.append(fn_ackley_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1]))
tuplaParticles = []
#initialize the particles
for i in range(1, PARTICLES):     
    particles.append(Particle())  
    #print(i)
    if(math.fabs(particles[i].xFitness) < bestAbsFitness):
        bestAbsFitness = math.fabs(particles[i].xFitness)
        V_GLOBAL_BEST[0] = particles[i].vectorpBest[0]
        V_GLOBAL_BEST[1] = particles[i].vectorpBest[1]
        V_GLOBAL_BEST[2] = fn_ackley_function(V_GLOBAL_BEST[0], V_GLOBAL_BEST[1])
 #   tuplaParticles.append([particles[i].vectorX[0], particles[i].vectorX[1]])
#part = (tuplaParticles)

fig = plt.figure()
x = []
y = []
for i in range(PARTICLES):
    x.append(particles[i].vectorX[0])
    y.append(particles[i].vectorX[1])

ax = fig.add_subplot(111)
ax.grid(True, linestyle = '-', color = '0.75')
ax.set_xlim([-3, 3])
ax.set_ylim([-3, 3])

scat = plt.scatter(x, y, c=x)
scat.set_alpha(0.8)

anim = animation.FuncAnimation(fig, _update_plot, fargs = (fig, scat),interval = 100)

plt.show()