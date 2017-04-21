from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors, ticker, cm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

data = {"om_dr":[], "freq":[], "apt": []} 

for line in open("frequency.txt"):
#	print(line[0:17]+"end"
    data["om_dr"].append(float(line[0:17]))
    data["freq"].append(float(line[18:40]))
    data["apt"].append(float(line[41:]))

#get the x and y axis data

m = len(np.unique(data["om_dr"]))
n = len(np.unique(data["freq"]))

dimension = (m,n)
print(dimension)

#reararnge the array
mo_dr = np.array(data["om_dr"]).reshape(dimension)
freq = np.array(data["freq"]).reshape(dimension)
apt = np.array(data["apt"]).reshape(dimension)

# grid the plot area
gs = gridspec.GridSpec(3, 2)

ax = plt.subplot(gs[0, :])
plt.contourf(mo_dr, freq, apt, locator=ticker.LogLocator(), cmap=cm.Blues_r)
plt.colorbar()
ax.set_title('whole')



ax0 = plt.subplot(gs[1,0])
plt.plot(freq[0],apt[0])
ax0.set_title('$\Omega_{dr}=0$')

ax1 = plt.subplot(gs[1,1])
plt.plot(freq[3],apt[3])
ax1.set_title('$\Omega_{dr}=0.06$')

ax2 = plt.subplot(gs[2,0])
plt.plot(freq[6],apt[2])
ax2.set_title('$\Omega_{dr}=0.12$')

ax3 = plt.subplot(gs[2,1])
plt.plot(freq[10],apt[10])
ax3.set_title('$\Omega_{dr}=0.2$')


plt.show()






