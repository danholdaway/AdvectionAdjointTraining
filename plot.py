import numpy as np
import matplotlib
import matplotlib.pyplot as plt

file1 = 'q_init.txt'
file2 = 'q_final.txt'

plotfile = 'q.png'

fs1 = [line.rstrip('\n') for line in open(file1)]
f1 = [float(x) for x in fs1]

fs2 = [line.rstrip('\n') for line in open(file2)]
f2 = [float(x) for x in fs2]

nx = len(fs1)
x = np.linspace(0,2*np.pi,nx)

plt.plot(x,f1,'b')
plt.plot(x,f2,'r--')
plt.savefig(plotfile)
plt.show()
