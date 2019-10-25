#  © Shahram Talei @ 2019 The University of Alabama
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
	coordinatesDM = ad[("Halo","Coordinates")]
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
	print("separations in contaminations:")
	print(r)
	# should be Rv instead of 0.17 but I have to check their unit and convert them.
	rContamination=r[r<(Rv)]
	count=len(rContamination)
	print("Halo %d at:" % (h.id))
	#print(center)
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
    def __init__(self,hid,R,M,a,b,c):
        self.hid=hid
        self.R=R
        self.M=M
        self.b_a=b/a
        self.c_a=c/a
        self.T=(a**2.-b**2.)/(a**2.+b**2.)
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
                    if(self.IsInside(point,Rout)):# and not(self.IsInside(point,Rin))):
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
    coordinatesDM = ad[("Halo","Coordinates")]
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
    coordsDM=np.array(coordinatesDM.v)
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
    Pmass=ad[("Halo","Mass")].in_units('Msun')
    print("Individual particle mass: %g"%Pmass[0])
    IDsDM = ad[("Halo","ParticleIDs")]
    print(len(IDsDM))
    hid =h.id# halo.quantities['particle_identifier']
    #halos = HaloFinder(ds, ptype=Halo, dm_only=False, total_mass=None)
    #ind = halos[0]["particle_index"] # list of particles IDs in this halo
    #print(len(ind))
    # REMOVE
    #Rvir=10
	# Rem
    bins=2.
    iteLim=3
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
            if (i==0):
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
        FullShape.append(shape(h.hid,h.R,h.M,axis[0],axis[1],axis[2]))
        #a[i]=axis[0] # use max instead
        #b[i]=axis[1]
        #c[i]=axis[2] # use min instead
        #b_a[i]=b[i]/a[i]
        #c_a[i]=c[i]/a[i]
        #use the shape object instead
    #now to save these numbers somewhere
    # save halo id, center, bins, and shapes,

    #and let's plot these




#how to run: python ShapeAnalysis.py snapshot_file halo_catalog check_contamination extract_shape
#example: $python ShapeAnalysis.py snap_264 halos_0.0.bin 1 1

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("snap", type=str)
	parser.add_argument("halo",type=str)
	parser.add_argument("contamination", type=int)
	parser.add_argument("extractshape", type=int)
	args = parser.parse_args()
	ds = yt.load(args.snap)#, unit_base=unit_base1)#,unit_system='galactic')
	if ds is None:
		print ("Error, sorry, I couldn't read the snapshot.!")
		sys.exit(1)
	print("Length unit: ", ds.length_unit.in_units('Mpc'))
	print("Time unit: ", ds.time_unit.in_units('Gyr'))
	print("Mass unit: ", ds.mass_unit.in_units('Msun'))
	print("Velocity unit: ", ds.velocity_unit.in_units('km/s'))
	#dh=yt.load(args.halo)
	dh=np.genfromtxt('halos_0.0.ascii', skip_header=18)
	if dh is None:
		print ("Error, sorry, I couldn't read the halo binary file.!")
		sys.exit(1)
	#halos = dh.all_data()
	#hc = HaloCatalog(data_ds=ds, finder_method='hop', finder_kwargs={ "dm_only": False , "ptype":"Halo"} )
	#hc = HaloCatalog(data_ds=ds, halos_ds=dh)#, finder_method='rockstar')#, finder_kwargs={ "dm_only": False })#, "dm_type":1} )
	#hc.add_callback("virial_quantities", ["radius"], profile_storage = "virial_quantities_profiles")
	#hc.add_recipe("calculate_virial_quantities", ["radius"])
	#hc.create(save_halos=True)
	#ad0 = hc.halos_ds.all_data()
	#masses = ad0['particle_mass'][:]#.in_units('Msun')
	#Rvir = ad0['virial_radius'][:]
	#Rvir=Rvir.in_units('Mpc')#*=3.24078e-25#=0.2
	#R2=halos["halos", "virial_radius"]
	#Rvir=0.25#175
	#print("Mvir=%g-Rvir:%g"%(masses,Rvir))
	Idh=np.array(dh[:,0])
	#CountAll= len(id)
	p=1000
	pnumh=np.array(dh[:,1])
	Mvirh=np.array(dh[:,2])
	Rvirh=np.array(dh[:,4])# in kpc
	xhalo=np.array(dh[:,8])
	yhalo=np.array(dh[:,9])
	zhalo=np.array(dh[:,10])
	#Posh=(xhalo,yhalo,zhalo)
	#print(Posh)

	print("Mvir-Rvir:")
	print(Mvirh)
	print(Rvirh)
	print(xhalo)
	print(pnumh)
	hid=Idh[pnumh>p]
	hM=Mvirh[pnumh>p]
	hx=xhalo[pnumh>p]
	hy=yhalo[pnumh>p]
	hz=zhalo[pnumh>p]
	hR=Rvirh[pnumh>p]
	hR/=1000 # convert to Mpc
	#print(hR)
	hP=pnumh[pnumh>p]
	print(hR)
	center=np.zeros(3)
	center[0]=hx
	center[1]=hy
	center[2]=hz
	print("# halos valid for the shape anlysis:%d"%len(hM))
	h=halo(center,hR,hM,hP,hid) # be careful later if you have more than one halo
	# check for contamination
	if(args.contamination == 1):
		print("checking for contamination")
		if(not(IsContaminated(h,ds))):
			print("halo is not contaminated")
		#add_callback("IsContaminated", _IsContaminated)
		#hc.add_callback("IsContaminated", ds, Rvir)
	#hc.create(save_halos=True)
	#
	if(args.extractshape == 1):
		print("Let's extract the shape!")
		GetShape(h,ds)
	#add_callback("GetShape", _GetShape)
	#hc.add_callback("GetShape",ds,Rvir)
	#hc.create(save_halos=True)
        #ad = ds.all_data()
    	#coordinatesDM = ad[("Halo","Coordinates")]
	#ellipsoid=a,b,c and rotation
	#add_callback(get_shape(coords,elipsoid,Rin,Rout))
	#ellipsoid=a,b,c and rotation
