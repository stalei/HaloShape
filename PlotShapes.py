#  © Shahram Talei @ 2020 The University of Alabama - All rights reserved.
#you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
import csv
import matplotlib.pyplot as plt
import numpy as np

#files
fG='/media/shahram/SD/19880/g3z19880-976.0.csv'
fC='/media/shahram/SD/19880/c3z19880-932.0.csv'
hID=19880

#open and plot

with open(fG, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    next(reader)
    next(reader)
    next(reader)
    size=8 # I have to add a better way to get this number
    print("Number of bins:%d"%size)
    R=[0]*size
    b_a=[0]*size
    c_a=[0]*size
    AJ=[0]*size
    AA0=[0.0]*size
    i=0
    for row in reader:
        R[i]=np.around(float(row[1]),3)
        b_a[i]=np.around(float(row[2]),3)
        AJ[i]=np.around(float(row[5]),3)
        AA0[i]=np.around(float(row[6]),3)
        i+=1
        print(row[1])
    fig1 = plt.figure(figsize=plt.figaspect(1./3.))
    ax11 = fig1.add_subplot(131)
    ax11.set_xlabel('R (Mpc)')
    ax11.set_ylabel('Axis Ratio (Mpc)')
    ax11.set_ylim(0,1)
    ax11.plot(R,b_a,'r',linestyle='-',label="b/a")#str(h.id))
    ax12 = fig1.add_subplot(132)
    ax12.set_xlabel('R (Mpc)')
    ax12.set_ylabel('Shape-Angular Momentum')
    ax12.set_ylim(0,1)
    ax12.plot(R,AJ,'b',linestyle='-.',label="NoSub")#+str(h.id))
    ax13 = fig1.add_subplot(133)
    ax13.set_xlabel('R (Mpc)')
    ax13.set_ylabel('Residual Orientation')
    ax13.set_ylim(0,1)
    ax13.plot(R,AA0,'b',linestyle='-.',label="NoSub")#+str(h.id))
    print(R)
    print(b_a)
    plt.show()
