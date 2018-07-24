import numpy as np
import matplotlib
import matplotlib.pyplot as plt

file1 = 'q_init_ad.txt'
file2 = 'q_final_ad.txt'
file3 = 'u_final_ad.txt'

plotfile = 'q.png'

fs1 = [line.rstrip('\n') for line in open(file1)]
f1 = [float(x) for x in fs1]

fs2 = [line.rstrip('\n') for line in open(file2)]
f2 = [float(x) for x in fs2]

fs3 = [line.rstrip('\n') for line in open(file3)]
f3 = [float(x) for x in fs3]

nx = len(fs1)
x = np.linspace(0,2*np.pi,nx)

plt.plot(x,f1,'b')
plt.plot(x,f2,'r--')
plt.plot(x,f3,'g--')
plt.savefig(plotfile)
plt.show()
