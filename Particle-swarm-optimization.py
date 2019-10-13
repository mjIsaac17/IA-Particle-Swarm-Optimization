import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

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

a = animation.FuncAnimation(fig, update, next, blit = False, interval = 10,
                               repeat = False)
plt.show()