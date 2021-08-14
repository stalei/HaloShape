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
from matplotlib.pyplot import text
import statistics
plt.rcParams["font.size"] =12
#files
#redshift 0
f1='/media/shahram/SD/Sample100Mpc/32251/g3z32251_377.0.csv'
f2='/media/shahram/SD/Sample100Mpc/32251/c3z32251_426.0.csv'

f3='/media/shahram/SD/Sample100Mpc/717/g3z717_2023.0.csv'
f4='/media/shahram/SD/Sample100Mpc/717/c3z717_2002.0.csv'

f5='/media/shahram/SD/Sample100Mpc/19880/g3z19880_976.0.csv'
f6='/media/shahram/SD/Sample100Mpc/19880/c3z19880_932.0.csv'
#redshift 3
f7='/media/shahram/SD/Sample100Mpc/717/z3/2888.0-G717-z3-2.csv'
f8='/media/shahram/SD/Sample100Mpc/717/z3/2899.0-C717z3.csv'

f9='/media/shahram/SD/Sample100Mpc/19880/z3/2658.0-G19880z3.csv'
f10='/media/shahram/SD/Sample100Mpc/19880/z3/2679.0-C19880z3.csv'

f11='/media/shahram/SD/Sample100Mpc/32251/z3/873.0-G32251z3.csv'
f12='/media/shahram/SD/Sample100Mpc/32251/z3/852.0-C32251z3.csv'

f13='/media/shahram/SD/Sample100Mpc/32251/High5M/9678.0.G.csv'

f15='/media/shahram/SD/Sample100Mpc/32251/High5M/z1.5/21162.0-G.csv'
f16='/media/shahram/SD/Sample100Mpc/32251/High5M/z1.5/22045.0-C.csv'

f17='/media/shahram/SD/Sample100Mpc/32251/High5M/z0.62/14003.0-G-30bin.csv'
f18='/media/shahram/SD/Sample100Mpc/32251/High5M/z0.62/15486.0-C-30bin.csv'

f19='/media/shahram/SD/Sample100Mpc/m12b/z0/6520-m12bGz0 (copy).csv'
f20='/media/shahram/SD/Sample100Mpc/m12b/z0/4743-m12bCz0 (copy).csv'

f21='/media/shahram/SD/Sample100Mpc/32251/High5M/z0.18/10096-32251Gz0.18.csv'
f22='/media/shahram/SD/Sample100Mpc/32251/High5M/z0.18/11607-32251Cz0.18.csv'

f23='/media/shahram/SD/Sample100Mpc/32251/High5M/z0/9678-32251Gz0.csv'
f24='/media/shahram/SD/Sample100Mpc/32251/High5M/z0/9575-32251Cz0.csv'

f25='/media/shahram/SD/Sample100Mpc/m12b/z0.3/7154-m12bGz0.3.csv'
f26='/media/shahram/SD/Sample100Mpc/m12b/z0.3/5526-m12bCz0.3.csv'

h_title='m12b,z=0'

#open and plot

with open(f20, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    next(reader)
    next(reader)
    next(reader)
    size=29 # I have to add a better way to get this number
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
        c_a[i]=np.around(float(row[3]),3)
        b_aNoSub[i]=np.around(float(row[5]),3)
        c_aNoSub[i]=np.around(float(row[6]),3)
        AJ[i]=np.around(float(row[8]),3)
        BJ[i]=np.around(float(row[10]),3)
        CJ[i]=np.around(float(row[12]),3)
        AA0[i]=np.around(float(row[9]),3)
        BB0[i]=np.around(float(row[11]),3)
        CC0[i]=np.around(float(row[13]),3)
        i+=1
        print(row[1])
    fig1 = plt.figure(1,figsize=plt.figaspect(1./3.))
    ax11 = fig1.add_subplot(131)
    ax11.set_xlabel('R (kpc)')
    ax11.set_ylabel('Axis Ratio')
    ax11.set_ylim(0,1.2)
    ax11.plot(R,b_a,'r',linestyle='-',label="b/a")#str(h.id))
    ax11.plot(R,c_a,'b',linestyle='-',label="c/a")#str(h.id))
    #ax11.plot(R,b_aNoSub,'r',linestyle='-.',label="b/a-No Subhalo")#str(h.id))
    #ax11.plot(R,c_aNoSub,'b',linestyle='-.',label="c/a-No Subhalo")#str(h.id))
    linestyles = ['-', ":"]
    dummy_lines = []
    for b_idx in range(0,2):# enumerate(np.unique(df["b"])):
        dummy_lines.append(ax11.plot([],[], c="black", ls = linestyles[b_idx])[0])
    legend2 = ax11.legend([dummy_lines[i] for i in [0,1]], ["CoSANG", "DMO"], loc=2)
    ax11.add_artist(legend2)
    ax11.legend(loc=1)
    #ax11.add_artist(legend2)
    #text(5.,1.15 , " ̶  CoSANG", rotation=0, verticalalignment='center')
    #text(5.,1.09 , ".. DMO", rotation=0, verticalalignment='center')
    ax12 = fig1.add_subplot(132)
    ax12.set_xlabel('R (kpc)')
    ax12.set_ylabel('Shape-Angular Momentum')
    ax12.set_ylim(0,1.2)
    ax12.plot(R,AJ,'b',linestyle='-',label="$\hat{A}.\hat{L}$")#+str(h.id))
    ax12.plot(R,BJ,'r',linestyle='-',label="$\hat{B}.\hat{L}$")#+str(h.id))
    ax12.plot(R,CJ,'g',linestyle='-',label="$\hat{C}.\hat{L}$")#+str(h.id))
    legend3 = ax12.legend([dummy_lines[i] for i in [0,1]], ["CoSANG", "DMO"], loc=1,ncol=2, bbox_to_anchor=(0.65,0.916))
    ax12.add_artist(legend3)
    ax12.legend(loc=2,ncol=3)
    ax13 = fig1.add_subplot(133)
    ax13.set_xlabel('R (kpc)')
    ax13.set_ylabel('Residual Orientation')
    ax13.set_ylim(0,1.2)
    ax13.plot(R,AA0,'b',linestyle='-',label="$\hat{A}.\hat{A}_0$")#+str(h.id))
    ax13.plot(R,BB0,'r',linestyle='-',label="$\hat{B}.\hat{B}_0$")#+str(h.id))
    ax13.plot(R,CC0,'g',linestyle='-',label="$\hat{C}.\hat{C}_0$")#+str(h.id))
    legend4 = ax13.legend([dummy_lines[i] for i in [0,1]], ["CoSANG", "DMO"], loc=4)
    ax13.add_artist(legend4)
    ax13.legend(loc=3)
    print(R)
    print(b_a)
    fig2=plt.figure(2,figsize=plt.figaspect(1.))
    ax21 = fig2.add_subplot(111)
    ax21.title.set_text(h_title)
    ax21.set_xlabel('b/a')
    ax21.set_ylabel('c/a')
    ax21.scatter(b_a,c_a,marker='o',color='b',label='CoSANG')
    #ax21.plot(b_a,b_a,'r',linestyle='-',label="Prolate")#str(h.id))
    text(0.66,0.696 , "Prolate", rotation=45, verticalalignment='center')
    #text(np.max(b_a)+0.2,np.max(c_a)-0.2 , "Oblate", rotation=90, verticalalignment='center')
    #text(0.63,0.78 , "Sphere", rotation=0, verticalalignment='center')
    b_a_ave=statistics.mean(b_a)
    c_a_ave=statistics.mean(c_a)
    b_a_std=np.std(b_a)
    c_a_std=np.std(c_a)
    ax21.errorbar(b_a_ave, c_a_ave, b_a_std,c_a_std, linestyle='None', marker='x',label='CoSANG')
    #plt.show()
####
# 2nd plot

with open(f19, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    next(reader)
    next(reader)
    next(reader)
    #size=15 # I have to add a better way to get this number
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
        c_a[i]=np.around(float(row[3]),3)
        b_aNoSub[i]=np.around(float(row[5]),3)
        c_aNoSub[i]=np.around(float(row[6]),3)
        AJ[i]=np.around(float(row[8]),3)
        BJ[i]=np.around(float(row[10]),3)
        CJ[i]=np.around(float(row[12]),3)
        AA0[i]=np.around(float(row[9]),3)
        BB0[i]=np.around(float(row[11]),3)
        CC0[i]=np.around(float(row[13]),3)
        i+=1
        print(row[1])
    fig1 = plt.figure(1,figsize=plt.figaspect(1./3.))
    #ax11 = fig1.add_subplot(131)
    ax11.set_xlabel('R (kpc)')
    ax11.set_ylabel('Axis Ratio')
    ax11.set_ylim(0,1.2)
    ax11.plot(R,b_a,'r',linestyle=':',label="b/a")#str(h.id))
    ax11.plot(R,c_a,'b',linestyle=':',label="c/a")#str(h.id))
    ax11.title.set_text(h_title)
    #ax11.plot(R,b_aNoSub,'r',linestyle='-.',label="b/a-No Subhalo")#str(h.id))
    #ax11.plot(R,c_aNoSub,'b',linestyle='-.',label="c/a-No Subhalo")#str(h.id))
    #ax11.legend(loc=1)
    #ax12 = fig1.add_subplot(132)
    #ax12.set_xlabel('R (kpc)')
    #ax12.set_ylabel('Shape-Angular Momentum')
    ax12.set_ylim(0,1.2)
    ax12.plot(R,AJ,'b',linestyle=':',label="$\hat{A}.\hat{L}$")#+str(h.id))
    ax12.plot(R,BJ,'r',linestyle=':',label="$\hat{B}.\hat{L}$")#+str(h.id))
    ax12.plot(R,CJ,'g',linestyle=':',label="$\hat{C}.\hat{L}$")#+str(h.id))
    ax12.title.set_text(h_title)
    #ax12.legend(loc=2,ncol=3)
    #ax13 = fig1.add_subplot(133)
    #ax13.set_xlabel('R (kpc)')
    #ax13.set_ylabel('Residual Orientation')
    ax13.set_ylim(0,1.2)
    ax13.plot(R,AA0,'b',linestyle=':',label="$\hat{A}.\hat{A}_0$")#+str(h.id))
    ax13.plot(R,BB0,'r',linestyle=':',label="$\hat{B}.\hat{B}_0$")#+str(h.id))
    ax13.plot(R,CC0,'g',linestyle=':',label="$\hat{C}.\hat{C}_0$")#+str(h.id))
    ax13.title.set_text(h_title)
    #ax13.legend(loc=3)
    print(R)
    print(b_a)
    fig2=plt.figure(2,figsize=plt.figaspect(1.))
    #ax21 = fig2.add_subplot(111)
    ax21.set_xlabel('b/a')
    ax21.set_ylabel('c/a')
    ax21.scatter(b_a,c_a,marker='s',color='r',label='DMO')
    ax21.plot(b_a,b_a,'k',linestyle='-')#str(h.id))
    #text(0.4,0.5 , "Prolate", rotation=45, verticalalignment='center')
    #text(0.77,0.4 , "Oblate", rotation=90, verticalalignment='center')
    text(np.max(b_a)+0.1,np.max(c_a)-0.3 , "Oblate", rotation=90, verticalalignment='center')
    #text(0.63,0.78 , "Sphere", rotation=0, verticalalignment='center')
    b_a_ave=statistics.mean(b_a)
    c_a_ave=statistics.mean(c_a)
    b_a_std=np.std(b_a)
    c_a_std=np.std(c_a)
    ax21.errorbar(b_a_ave, c_a_ave, b_a_std,c_a_std, linestyle='None', marker='x',label='DMO')
    ax21.legend(loc=2)
    plt.show()
