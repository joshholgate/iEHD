import mayavi.mlab as mlab
import numpy as np
import sys
#import matplotlib.pyplot as plt
from scipy import interpolate

data = np.loadtxt("interface.dat", usecols=(0,1))
data=data[np.argsort(data[:,0])]

r = interpolate.interp1d(data[:,0],data[:,1],)
zmin = data[0,0]
zmax = data[-1,0]

#z = np.arange(zmin, zmax, 0.01)
#plt.plot(data[:,0],data[:,1], 'ro')
#plt.plot(z,r(z),'b-')
#plt.show()

dz, dtheta = 0.01, np.pi/250.0
[z,theta] = np.mgrid[zmin:zmax:dz,0:2*np.pi+dtheta*1.5:dtheta]
x = r(z)*np.cos(theta)
y = r(z)*np.sin(theta)

# View it.
mlab.figure(bgcolor=(1, 1, 1))
s = mlab.mesh(x, y, z, color=(0.5,0.5,0.5))
mlab.view(azimuth=0, elevation=90, distance=8, roll=0)
mlab.savefig("t"+str(sys.argv[1])+"_3d_drop.png",size=(800,636))
