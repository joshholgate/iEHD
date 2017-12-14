import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import interpolate
from matplotlib.mlab import griddata

data = np.loadtxt("t"+str(sys.argv[1])+"_gridvalues.dat", usecols=(0,1,5))

x = data[:,1]
y = data[:,0]
LS = data[:,2]
xi = np.linspace(0.0, 6.0, 200)
yi = np.linspace(0.0, 2.0, 200)
# grid the data.
LSi = griddata(x, y, LS, xi, yi, interp='linear')

CS = plt.contour(xi, yi, LSi, levels=[0])

p = CS.collections[0].get_paths()[0]
v = p.vertices

np.savetxt("interface.dat",v)
