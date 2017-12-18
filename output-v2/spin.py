import numpy as np
#import matplotlib.pyplot as plt
#from scipy import interpolate
from scipy import ndimage

open("crit_points.dat", "w").close()

fissility = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

for i in range(0, len(fissility)):
  data = np.loadtxt("x"+str(fissility[i])+".dat", usecols=(1,3))
  omega_array = ndimage.filters.gaussian_filter(data[:,1],3000)
  crit_point = np.argmax(omega_array)
  with open("crit_points.dat", "a") as myfile:
    myfile.write(str(fissility[i])+"\t"+str(data[crit_point,0])+"\t"+str(data[-1,0])+"\n")

#plt.plot(data[:,0],data[:,1], 'r-')
#plt.plot(l,omega(l),'b-')
#plt.show()
with open("crit_points.dat", "a") as myfile:
  myfile.write("1\t0\t0")
