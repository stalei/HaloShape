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
plt.rcParams["font.size"] =12
#files
f1='/media/shahram/SD/19880/g3z19880-976.0.csv'
f2='/media/shahram/SD/19880/c3z19880-932.0.csv'
f3='/media/shahram/SD/717/c3z19880-932.0.csv'
f5='/media/shahram/SD/717/c3z19880-932.0.csv'
f6='/media/shahram/SD/32251/c3z19880-932.0.csv'
f7='/media/shahram/SD/32251/c3z19880-932.0.csv'

#open and plot

with open(f1, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    next(reader)
    next(reader)
    next(reader)
    size=8 # I have to add a better way to get this number
    print("Number of bins:%d"%size)
    #declare variables
    R=[0]*size #for all plots
    #sub 131
    b_a=[0]*size
    b_aNoSub=[0]*size
    c_a=[0]*size
    c_aNoSub=[0]*size
    #sub 132
    AJ=[0]*size
    BJ=[0]*size
    CJ=[0]*size
    #sub 133
    AA0=[0.0]*size
    BB0=[0.0]*size
    CC0=[0.0]*size

    i=0
    for row in reader:
        R[i]=np.around(float(row[1]),3)*1000
        b_a[i]=np.around(float(row[2]),3)
        AJ[i]=np.around(float(row[5]),3)
        AA0[i]=np.around(float(row[6]),3)
        i+=1
        print(row[1])
    fig1 = plt.figure(figsize=plt.figaspect(1./3.))
    ax11 = fig1.add_subplot(131)
    ax11.set_xlabel('R (kpc)')
    ax11.set_ylabel('Axis Ratio')
    ax11.set_ylim(0,1)
    ax11.plot(R,b_a,'r',linestyle='-',label="b/a")#str(h.id))
    ax11.plot(R,c_a,'b',linestyle='-',label="c/a")#str(h.id))
    ax11.plot(R,b_aNoSub,'r',linestyle='-.',label="b/a-No Subhalo")#str(h.id))
    ax11.plot(R,c_aNoSub,'b',linestyle='-.',label="c/a-No Subhalo")#str(h.id))
    ax11.legend(loc=4)
    ax12 = fig1.add_subplot(132)
    ax12.set_xlabel('R (kpc)')
    ax12.set_ylabel('Shape-Angular Momentum')
    ax12.set_ylim(0,1)
    ax12.plot(R,AJ,'b',linestyle='-',label="$\hat{A}.\hat{L}$")#+str(h.id))
    ax12.plot(R,BJ,'r',linestyle='-.',label="$\hat{B}.\hat{L}$")#+str(h.id))
    ax12.plot(R,CJ,'g',linestyle=':',label="$\hat{C}.\hat{L}$")#+str(h.id))
    ax12.legend(loc=1)
    ax13 = fig1.add_subplot(133)
    ax13.set_xlabel('R (kpc)')
    ax13.set_ylabel('Residual Orientation')
    ax13.set_ylim(0,1)
    ax13.plot(R,AA0,'b',linestyle='-',label="$\hat{A}.\hat{A}_0$")#+str(h.id))
    ax13.plot(R,BB0,'r',linestyle='-.',label="$\hat{B}.\hat{B}_0$")#+str(h.id))
    ax13.plot(R,CC0,'g',linestyle=':',label="$\hat{C}.\hat{C}_0$")#+str(h.id))
    ax13.legend(loc=1)
    print(R)
    print(b_a)
    plt.show()
