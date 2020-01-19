# Double check the moment tensor

#  © Shahram Talei @ 2019 The University of Alabama - All rights reserved.
#you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
# the main differnce in V2 is the COM instead of the center coordinates from rockstar

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import linalg as LA

x=np.genfromtxt("x.out",delimiter=',',comments='#')
y=np.genfromtxt("y.out",delimiter=',',comments='#')
z=np.genfromtxt("z.out",delimiter=',',comments='#')
fig0 = plt.figure(0,figsize=plt.figaspect(1))
ax0 = fig0.add_subplot(111, projection='3d')
ax0.scatter(x,y,z,c='black',alpha=0.9,marker='.',s=1)
coords=np.zeros((len(x),3))
coords[:,0]=x
coords[:,1]=y
coords[:,2]=z
shape=[[0,0,0],[0,0,0],[0,0,0]]

for i in range(0,3):
    for j in range(0,3):
        c=s=0
        print("s[%d,%d]"%(i,j))
        for point in coords:
            print("point before if:")
            print(point)
            print("pp/s")
            print(point[i]*point[j])
            s+=point[i]*point[j]
            print(s)
            c+=1
        #Mtot=c*m
        if c>0:
            shape[i][j]=s/c
print(shape)
s2=np.array(shape)
sfull=s2.T+s2
for i in range(0,3):
    sfull[i][i]/=2.
print(sfull-shape)
ax,ori=LA.eig(shape)
print(ax)
print(ori)

ax2,ori2=LA.eig(sfull)
print(ax2)
print(ori2)



plt.show()



# a=[[0.00495025,1.77223319313493e-020,-1.92592994438724e-035],[1.77223319313493e-020,0.00297015,1.04366941327238e-035],[-1.92592994438724e-035,1.04366941327238e-035,0.0118806]]
#>>> from scipy import linalg as LA
#>>> val,vec=LA.eig(a)
#>>> val
#array([0.00495025+0.j, 0.00297015+0.j, 0.0118806 +0.j])
#>>> vec
#array([[ 1.00000000e+00, -8.95022066e-18,  2.77897934e-33],
#       [ 0.00000000e+00,  1.00000000e+00,  0.00000000e+00],
#       [ 0.00000000e+00,  0.00000000e+00, -1.00000000e+00]])
