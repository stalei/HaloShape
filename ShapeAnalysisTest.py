#  © Shahram Talei @ 2019 The University of Alabama - All rights reserved.
#you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division
import yt
import numpy as np
from yt.analysis_modules.halo_finding.api import *
from yt.analysis_modules.halo_analysis.api import *
from os import environ
environ['CFLAGS'] = "-I"+np.get_include()

import pyximport; pyximport.install()
#import particle_ops
import argparse


import tempfile
import shutil
import os
import sys

from scipy.spatial.transform import Rotation as R
from numpy import linalg as LA
from operator import mul
from functools import reduce
import matplotlib.pyplot as plt
from yt.units import parsec, Msun, Mpc

#define an ellipsoid object in an arbitrary orientation (a,b,c,Alpha,Beta,Yamma)

unit_base0 = {'UnitLength_in_cm'         : 3.08568e+24,
             'UnitMass_in_g'            :   1.989e+43,
             'UnitVelocity_in_cm_per_s' :      100000}

unit_base1 = {'length': (1.0, 'Mpc'),
    'velocity': (1.0, 'km/s'),
    'mass': (1.0e10, 'Msun')}

def IsContaminated(h,ds):#snap):
	#halo
	#dds = halo.halo_catalog.data_ds
	#print(halo.quantities[""])
	center =h.pos# dds.arr([halo.quantities["particle_position_%s" % axis] \
	Rv=h.R
	#print(h.R)
	#for axis in "xyz"])
	#my_id = halo.quantities['particle_identifier']
	#center*=3.24078e-25
	#print("Halo %d" % (my_id))
	#paarticles
	count=0
	boundary=np.zeros((2,3))
	#ds = yt.load(args.snap)#,index_ptype="Halo")
	#dsBND = yt.load(args.snap,index_ptype="Bndry")
	ad = ds.all_data()
	#print(ad[("halo","ParPosX")])
	coordinatesDM = ad[("all","particle_position_x")]
	#massesDM=ad[("Halo","Masses")]
	coordinatesLowRes = ad[("Bndry","Coordinates")]
	DMCount=len(coordinatesDM)
	LowResCount=len(coordinatesLowRes)
	print("number of High resolution particles:%d"%DMCount)
	print("number of Low resolution particles:%d"%LowResCount)
	#print(len(coordinatesDM[:,0])), x for all particles
	for i in range(0,3):
		boundary[0,i]=np.min(coordinatesDM[:,i]) # (min max, x y z)
		boundary[1,i]=np.max(coordinatesDM[:,i])
	print("high resolution particles are in box:")
	print(boundary)
	#for i in range(0,LowResCount): # this for is super slow!!
	#	x=coordinatesLowRes[i,0]
	#	y=coordinatesLowRes[i,1]
	#	z=coordinatesLowRes[i,2]
	#	if(x>boundary[0,0] and x<boundary[1,0] and y>boundary[0,1] and y<boundary[1,1] and z>boundary[0,2] and z<boundary[1,2]):
	#		print("contamination at: %f,%f,%f"%x%y%z)
	#		count+=1
	xLR=coordinatesLowRes[:,0]
	yLR=coordinatesLowRes[:,1]
	zLR=coordinatesLowRes[:,2]

	xLR2=xLR[(xLR>boundary[0,0]) & (xLR<boundary[1,0])]
	yLR2=yLR[(xLR>boundary[0,0]) & (xLR<boundary[1,0])]
	zLR2=zLR[(xLR>boundary[0,0]) & (xLR<boundary[1,0])]

	xLR3=xLR2[(yLR2>boundary[0,1]) & (yLR2<boundary[1,1])]
	yLR3=yLR2[(yLR2>boundary[0,1]) & (yLR2<boundary[1,1])]
	zLR3=zLR2[(yLR2>boundary[0,1]) & (yLR2<boundary[1,1])]


	xLR4=xLR3[(zLR3>boundary[0,2]) & (zLR3<boundary[1,2])]
	yLR4=yLR3[(zLR3>boundary[0,2]) & (zLR3<boundary[1,2])]
	zLR4=zLR3[(zLR3>boundary[0,2]) & (zLR3<boundary[1,2])]

	#print(len(xLR7))
	#count=len(xLR7)
	r2=((xLR4.v-center[0])**2.+(yLR4.v-center[1])**2.+(zLR4.v-center[2])**2.)
	r=np.sqrt(r2)
	#print("separations in contaminations:")
	#print(r)
	# should be Rv instead of 0.17 but I have to check their unit and convert them.
	rContamination=r[r<(Rv)]
	count=len(rContamination)
	print("Halo %d at:" % (h.id))
	print(center)
	print("is contaminated with: %d low resolution particles at δr="%count)
	print(rContamination)
	return count

class halo:
    def __init__(self,pos,R,M,pnum,id):
        self.pos=pos
        self.R=R
        self.M=M
        self.pnum=pnum
        self.id=id
class shape:
    def __init__(self,hid,Rv,Mv,R,ax):
        ax.sort()
        self.hid=hid
        self.Rv=Rv
        self.Mv=Mv
        self.R=R
        self.a=ax[2]
        self.b=ax[1]
        self.c=ax[0]
        self.b_a=ax[1]/ax[2]
        self.c_a=ax[0]/ax[2]
        self.T=(ax[2]**2.-ax[1]**2.)/(ax[2]**2.-ax[0]**2.)
class ellipsoid:
    def __init__(self,axis,orientations):#put in tuple (self,p,angle)..p.x.
        self.a=axis[0]
        self.b=axis[1]
        self.c=axis[2]
        self.axis=axis
        self.orientations=orientations
        self.A=orientations[0]
        self.B=orientations[1]
        self.C=orientations[2]
    def IsInside(self, point,R):# put in tuple p ... p.x
        #r = R.from_euler('zyx', [self.A, self.B, self.C], degrees=True)
        #point=r.apply(point)
        #print("point in is inside:")
        #print(point)
        #print(R)
        t=0
        pNew=np.array([0.0,0.0,0.0])
        for i in range(0,3):
            for j in range(0,3):
                #print(point[j])
                #print(self.orientations[j,i])
                pNew[i]+=point[j]*self.orientations[j,i] # each eig vec is [:,i] comp
        #x=point[0]#-center[0]
        #y=point[1]#-center[1]
        #z=point[2]#-center[2]
        #print("point in is inside:")
        #print(pNew)
        #t=np.sum((pNew/self.axis)**2.)
        for i in range(0,3):
            t+=(pNew[i]/self.axis[i])**2.
        #t=(x/self.a)**2.+(y/self.b)**2.+(z/self.c)**2
        #print(t)
        if t<=1:
            return True
        else:
            return False
    def rotate(self, orientation):
        self.orientations=oriantation
    def volume(self):
        return self.axis[0]*self.axis[1]*self.axis[2]
    #    eulerrotation
    def ShapeTesnsor(self,coords,m, Rin,Rout):
        print("Started extracting the shape tensor")
        i=j=0
        c=0
        #N=2#len() use mass instead, you won't need to count!
        shape=[[0,0,0],[0,0,0],[0,0,0]]#np.array([[0,0,0,],[0,0,0,],[0,0,0]])
        s=0
        #points=coords[(coord[not IsInside()]) && ]
        #for i in range(0,3):
        #    p[i]=coords[:,i]
        #rEll=1#np.sqrt(c)Zemp 10
        #we can use IsInside to find both inside and outside
        for i in range(0,3):
            for j in range(0,3):
                c=s=0
                print("s[%d,%d]"%(i,j))
                for point in coords:
                    #print("point before if:")
                    #print(point)
                    if(self.IsInside(point,Rout)):# and (not(self.IsInside(point,Rin)))):
                        #print("point after if:")
                        #print(point)
                        #print("point:%g"%(m*point[i]*point[j]))
                        s+=point[i]*point[j]
                        c+=1
                #s[i,j]/=(N-1)*m
                Mtot=c*m
                if c>1:
                    shape[i][j]=s/(c-1)
                print("particle count for %g<R<%g=%d"%(Rin,Rout,c))
        return shape
def CompareEllipsoids(ell1,ell2):
    AreSame=False;
    error=1.0e-3
    c=0
    inP=np.array([0,0,0])
    dRatio=np.array([0,0])
    #dax=np.abs(np.divide(ell1.axis/ell2.axis))
    for i in range(1,3):
        if(ell1.axis[0] !=0):
            dRatio[i-1]=np.abs((ell1.axis[i]/ell1.axis[0])/(ell2.axis[i]/ell2.axis[0]))
        #dRatio[i-1]=0
    #dori=np.abs(np.divide(ell1.orientations/ell2.orientations))
    d=np.sum(dRatio[(dRatio<1+error) & (dRatio>1-error)])
    if(d==0):
        for i in range(0,3):
            for j in range(0,3):
                inP[i]+=ell1.orientations[j,i]*ell2.orientations[j,i]
            if(inP[i]<((ell1.axis[i])**2.+error) and inP[i]>((ell1.axis[i])**2.-error)):
                c+=1
        if(c==3):
            AreSame=True
    return AreSame;
def GetShape(h,ds):
    FullShape=[]
    ad = ds.all_data()
    xp=ad[("all","particle_position_x")]
    yp=ad[("all","particle_position_y")]
    zp=ad[("all","particle_position_z")]
    coordsArray=np.array((xp,yp,zp))
    coordinatesDM =coordsArray.T# we need to reshape
    #print(coordinatesDM)
    #dds = halo.halo_catalog.data_ds
    #center = dds.arr([halo.quantities["particle_position_%s" % axis] \
    #for axis in "xyz"])
    c2=r2=0
    #center*=3.24078e-25
    #center=center.in_units('Mpc')
    center =h.pos# dds.arr([halo.quantities["particle_position_%s" % axis] \
    Rvir=h.R
    print("halo center (Mpc):")
    print(center)
    coordsDM=np.array(coordinatesDM)
    for i in range(0,3):
        print(coordsDM[:,i])
        print(center[i])
        coordsDM[:,i]-=center[i]
        print(coordsDM[:,i])
        #print(i)
    #print(coordinatesDM)
    for i in range(0,3):
        r2+=(coordsDM[:,i])**2.
        #c2+=(center[i])**2.
    r=np.sqrt(r2)
    print(r)
    #coords=coordinatesDM[np.abs(r-np.sqrt(c2)<Rvir)]
    coords=coordsDM[r<=Rvir]
    print("# of virialized particles:%d"%len(coords))
    print(coords)
    #coords[:,0]-=center[0]
    #coords[:,1]-=center[1]
    #coords[:,2]-=center[2]
    #print(coords[0,:])
    Pmass=ad[("all","particle_mass")].in_units('Msun')
    print("Individual particle mass: %g"%Pmass[0])
    #IDsDM = ad[("Halo","ParticleIDs")]
    #print(len(IDsDM))
    hid =h.id# halo.quantities['particle_identifier']
    #halos = HaloFinder(ds, ptype=Halo, dm_only=False, total_mass=None)
    #ind = halos[0]["particle_index"] # list of particles IDs in this halo
    #print(len(ind))
    # REMOVE
    #Rvir=10
	# Rem
    bins=8
    iteLim=4
    #convLim=5 nor need, we just compare two
    #Rbins=np.logspace(0,Rvir,bins)#(Rvir/bins,Rvir,bins)
    Rbins=np.linspace(0,Rvir,bins+1)
    print("Rbins:")
    print(Rbins)
    axis=np.array([Rbins[1],Rbins[1],Rbins[1]])
    orientation=np.identity(3)#zeros((3,3))#.array([0,0,0],[0,0,0],[0,0,0])
    a=[0]*len(Rbins)
    b=[0]*len(Rbins)
    c=[0]*len(Rbins)
    b_a=[0]*len(Rbins)
    c_a=[0]*len(Rbins)
    v0=1#ell.volume
    for i in range(0,len(Rbins)-1):
        iteration=0
        convergence=False
        conv=0
        while(not(convergence) and iteration<iteLim):
            s=np.array([[1,0,0,],[0,1,0,],[0,0,1]])
            print("halo %d -bin %d -iteration %d"% (hid,i+1, iteration+1))
            ell=ellipsoid(axis,orientation)
            s=ell.ShapeTesnsor(coords,Pmass[0],Rbins[i],Rbins[i+1])
            if (iteration==0):
                v0=ell.volume()
            print("shape tensor:")
            print(np.array(s))
            axisNew,orientationNew=LA.eig(s)
            print("eigenvalues & eigenvectors:")
            print(axisNew)
            print(orientationNew)
            #axisNew=np.array([1,1,1]) #1 just for test
            #orientationNew=np.array([0,0,0],[0,0,0],[0,0,0])#orientationNew # 0 just for test
            #if() what if there was ne eignevalue?
            newell=ellipsoid(axisNew,orientationNew)
            if iteration>0:
                comp=CompareEllipsoids(ell,newell)
                if comp==True:
                    conv+=1
                    if conv>=1:
                        convergence=True
                else:
                    conv=0
                #if conv>=2:
                #    convergence=True#CompareEllipsoids(ell,newell)
                if(not convergence):
                    #iteration+=1  #  new_list = [x+1 for x in my_list]
                    vol=reduce(mul,axisNew)
                    norm=(v0/vol)**(1./3.)
                    axis=[ax*norm for ax in axisNew]
                    print("normalized axis:")
                    print(axis)
                    orientation=orientationNew
                if(iteration>=iteLim):
                    print("S didn't converge for halo %d"%hid)
                    convergence=True
            else:
                vol=reduce(mul,axisNew)
                print(vol)
                norm=(v0/vol)**(1./3.)
                axis=[ax*norm for ax in axisNew]
                orientation=orientationNew
                #if (i==0):
                #    v0=ell.volume()
            iteration+=1
        axis.sort()
        Rave=(Rbins[i]+Rbins[i+1])/2.
        FullShape.append(shape(h.id,h.R,h.M,Rave,axis))
        #a[i]=axis[0] # use max instead
        #b[i]=axis[1]
        #c[i]=axis[2] # use min instead
        #b_a[i]=b[i]/a[i]
        #c_a[i]=c[i]/a[i]
        #use the shape object instead
    #now to save these numbers somewhere
    # save halo id, center, bins, and shapes,
    print(FullShape[0].R)
    return FullShape

    #and let's plot these




#how to run: python ShapeAnalysis.py snapshot_file halo_catalog check_contamination extract_shape
#example: $python ShapeAnalysisTest.py snap_264 halos_0.0.bin 1 1

if __name__ == "__main__":
	p=1000
	UpperMass=1.0e13
	LowerMass=7.0e11

	center=np.zeros(3)
	center[0]=0
	center[1]=0
	center[2]=0
	n_particles = 30000
	h=halo(center,5.0e2,1.2e12,n_particles,0)
	ppx, ppy, ppz =1e2*np.random.normal(size=[3, n_particles])
	#print(np.shape(ppx))
	#ppx=range(0,10000,n_particles)
	#ppy=range(0,10000,n_particles)
	#ppz=range(0,10000,n_particles)
	ppm = np.ones(n_particles)
	data = {'particle_position_x': ppx,'particle_position_y': ppy,'particle_position_z': ppz,'particle_mass': ppm}
	bbox = 1.1*np.array([[min(ppx), max(ppx)], [min(ppy), max(ppy)], [min(ppz), max(ppz)]])
	#ds = yt.load_particles(data, length_unit=parsec, mass_unit=1e8*Msun, n_ref=256, bbox=bbox)
	ds = yt.load_particles(data, length_unit=Mpc, mass_unit=Msun, n_ref=256, bbox=bbox)
	#print(ds.field_list)
	print(ppx)

	fig2 = plt.figure(2)
	fig2.suptitle('Shapes')
	ax1 = fig2.add_subplot(221)
	ax2 = fig2.add_subplot(222)
	ax3 = fig2.add_subplot(223)
	ax4 = fig2.add_subplot(224)
	ax1.legend(loc=2)
	FinalHaloShape=[]
	print("Let's extract the shape!")
	HShape=GetShape(h,ds)
	ba=np.array(HShape[0].b_a)
	print(ba)
	pR=[]
	pb_a=[]
	pid=[]
	pc_a=[]
	pT=[]
	for s in HShape:
		#ax1.plot(s.R,s.b_a) or ax1.scatter(s.R,s.b_a,c='black', alpha=0.9, marker='.',s=15)
		#ax1.scatter( (s.R)*1000,s.b_a, alpha=0.9, marker='.',s=15)
		pR.append(s.R*1000)
		pb_a.append(s.b_a)
		pc_a.append(s.c_a)
		pid.append(s.hid)
		pT.append(s.T)
	ax1.title.set_text('b/a')
	ax2.title.set_text('c/a')
	ax3.title.set_text('T')
	ax4.title.set_text('c/a - b/a')
	#ax2.plot(pR,pb_a,pR,pb_a,'o',label=str(pid))
	ax1.plot(pR,pb_a,label=str(pid[0]))
	ax2.plot(pR,pc_a,label=str(pid[0]))
	ax3.plot(pR,pT,label=str(pid[0]))
	ax4.plot(pc_a,pb_a,'o',label=str(pid[0]))
	ax1.legend(loc=3)
	ax2.legend(loc=3)
	ax3.legend(loc=3)
	ax4.legend(loc=4)
	fig3=plt.figure(3)
	rp=np.sqrt(ppx**2.+ppy**2.+ppz**2.)
	plt.hist(rp,bins=100)
	plt.show()
