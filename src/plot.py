from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy

sol = numpy.loadtxt(open("solution.dat","rb"),delimiter=",")

x,y,z= sol[:,0], sol[:,1], sol[:,2]

fig = plt.figure()
ax = Axes3D(fig)

ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.2)
plt.show()
