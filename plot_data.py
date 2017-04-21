from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

data = {"om_dr":[], "freq":[], "apt": []} 
#data =[];
i = 1;

for line in open("frequency.txt"):
#	print(line[0:17]+"end"
    data["om_dr"].append(float(line[0:17]))
    data["freq"].append(float(line[18:40]))
    data["apt"].append(float(line[41:]))



x = data["om_dr"]
y = data["freq"]
z = data["apt"] 

plt.subplot(2, 1, 1)
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.2)

plt.xlabel('om_dr -- driving beam')
plt.ylabel('Frequency') 
#plt.zlabel('Amplitude')
plt.title('calculation _ fu')




plt.show()






