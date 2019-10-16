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


#define an ellipsoid object in an arbitrary orientation (a,b,c,Alpha,Beta,Yamma)

unit_base0 = {'UnitLength_in_cm'         : 3.08568e+24,
             'UnitMass_in_g'            :   1.989e+43,
             'UnitVelocity_in_cm_per_s' :      100000}

unit_base1 = {'length': (1.0, 'Mpc'),
    'velocity': (1.0, 'km/s'),
    'mass': (1.0e10, 'Msun')}

def _IsContaminated(halo,ds, Rv):#snap):
	#halo
	dds = halo.halo_catalog.data_ds
	#print(halo.quantities[""])
	center = dds.arr([halo.quantities["particle_position_%s" % axis] \
	for axis in "xyz"])
	my_id = halo.quantities['particle_identifier']
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
	r2=((xLR4-center[0])**2.+(yLR4-center[1])**2.+(zLR4-center[2])**2.)
	r=np.sqrt(r2)
	# should be Rv instead of 0.17 but I have to check their unit and convert them.
	rContamination=r[r<(Rv)]
	count=len(rContamination)
	print("Halo %d at:" % (my_id))
	print(center)
	print("is contaminated with: %d low resolution particles at δr="%count)
	print(rContamination)
	return count

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
        t=0
        pNew=np.array([0.0,0.0,0.0])
        for i in range(0,3):
            for j in range(0,3):
                #print(point[j])
                #print(self.orientations[j,i])
                pNew+=point[j].v*self.orientations[j,i] # each eig vec is [:,i] comp
        #x=point[0]#-center[0]
        #y=point[1]#-center[1]
        #z=point[2]#-center[2]
        #print("point in is inside:")
        #print(point)
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
        self.A=EulerAngles[0]
        self.B=EulerAngles[1]
        self.C=EulerAngles[2]
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
        rEll=1#np.sqrt(c)Zemp 10
        #we can use IsInside to find both inside and outside
        for i in range(0,3):
            for j in range(0,3):
                c=s=0
                print("s[%d,%d]"%(i,j))
                for point in coords:
                    if(self.IsInside(point,Rout)):# and not(self.IsInside(point,Rin))):
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
    dax=np.abs(np.divide(ell1.axis/ell2.axis))
    dori=np.abs(np.divide(ell1.orientations/ell2.orientations))
    sum[j]=ell1.orientations[:,i]*ell2.orientations[:,i]
    #uax=np.ones(3)
    #uori=np.ones(3,3)
    dax2=dax[dax<1.0001 & dax>0.9999]
    dori2=dori[dori<1.0001 & dori>0.9999]

    if(np.sum(dax2) ==0 and np.sum(dori2)==0):
        AreSame=True
    return AreSame;
def _GetShape(halo,ds,Rvir):
    ad = ds.all_data()
    coordinatesDM = ad[("Halo","Coordinates")]
    dds = halo.halo_catalog.data_ds
    center = dds.arr([halo.quantities["particle_position_%s" % axis] \
    for axis in "xyz"])
    c2=r2=0
    #center*=3.24078e-25
    center=center.in_units('Mpc')
    print("halo center (converted):")
    print(center)
    for i in range(0,3):
        r2+=(coordinatesDM[:,i]-center[i])**2.
        c2+=(center[i])**2.
    r=np.sqrt(r2)
    coords=coordinatesDM[np.abs(r-np.sqrt(c2)<Rvir)]
    print("# of virialized particles:%d"%len(coords))
    print(coords)
    for i in range(0,3):
        coords[:,i]-=center[i]
    #coords[:,0]-=center[0]
    #coords[:,1]-=center[1]
    #coords[:,2]-=center[2]
    print(coords[0,:])
    Pmass=ad[("Halo","Mass")].in_units('Msun')
    print("Individual particle mass: %g"%Pmass[0])
    IDsDM = ad[("Halo","ParticleIDs")]
    print(len(IDsDM))
    hid = halo.quantities['particle_identifier']
    #halos = HaloFinder(ds, ptype=Halo, dm_only=False, total_mass=None)
    #ind = halos[0]["particle_index"] # list of particles IDs in this halo
    #print(len(ind))
    # REMOVE
    #Rvir=10
	# Rem
    bins=2.
    iteLim=3
    #convLim=5 nor need, we just compare two
    Rbins=np.logspace(Rvir/bins,Rvir,bins)
    axis=np.array([Rbins[0],Rbins[0],Rbins[0]])
    orientation=np.identity(3)#zeros((3,3))#.array([0,0,0],[0,0,0],[0,0,0])
    a=[0]*len(Rbins)
    b=[0]*len(Rbins)
    c=[0]*len(Rbins)
    b_a=[0]*len(Rbins)
    c_a=[0]*len(Rbins)
    for i in range(0,len(Rbins)):
        iteration=0
        convergence=False
        conv=0
        while(not(convergence) and iteration<iteLim):
            s=np.array([[0,0,0,],[0,0,0,],[0,0,0]])
            print("halo %d -bin %d -iteration %d"% (hid,i+1, iteration+1))
            ell=ellipsoid(axis,orientation)
            s=ell.ShapeTesnsor(coords,Pmass[0],Rbins[i-1],Rbins[i])
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
            if i>0:
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
                    #iteration+=1
                    axis=axisNew
                    orientation=orientationNew
                if(iteration>=iteLim):
                    print("S didn't converge for halo %d"%hid)
                    convergence=True
            else:
                axis=axisNew
                orientation=orientationNew
                v0=ell.volume
            iteration+=1
        a[i]=axis[0]
        b[i]=axis[1]
        c[i]=axis[2]
        b_a[i]=b[i]/a[i]
        c_a[i]=c[i]/a[i]
    #now to save these numbers somewhere
    # save halo id, center, bins, and shapes,

    #and let's plot these




#python ShapeAnalysis.py snapshot_file check_contamination
#python ShapeAnalysis.py snap_264 1

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("snap", type=str)
	parser.add_argument("contamination", type=int)
	args = parser.parse_args()
	ds = yt.load(args.snap, unit_base=unit_base1,unit_system='galactic')
	if ds is None:
		print ("Error, sorry, I couldn't read the snapshot.!")
		sys.exit(1)
	print("Length unit: ", ds.length_unit.in_units('Mpc'))
	print("Time unit: ", ds.time_unit.in_units('Gyr'))
	print("Mass unit: ", ds.mass_unit.in_units('Msun'))
	print("Velocity unit: ", ds.velocity_unit.in_units('km/s'))
	hc = HaloCatalog(data_ds=ds, finder_method='hop', finder_kwargs={ "dm_only": False , "ptype":"Halo"} )
	#hc.add_callback("virial_quantities", ["radius"], profile_storage = "virial_quantities_profiles")
	#hc.add_recipe("calculate_virial_quantities", ["radius"])
	hc.create(save_halos=True)

	ad0 = hc.halos_ds.all_data()
	masses = ad0['particle_mass'][:].in_units('Msun')
	Rvir = ad0['virial_radius'][:]
	Rvir=Rvir.in_units('Mpc')#*=3.24078e-25#=0.2
	Rvir=0.175
	print("Mvir=%g-Rvir:%g"%(masses,Rvir))
	# check for contamination
	if(args.contamination == 1):
		print("checking for contamination")
		add_callback("IsContaminated", _IsContaminated)
		hc.add_callback("IsContaminated", ds, Rvir)
	#hc.create(save_halos=True)
	#
	#
	print("Let's extract the shape!")
	add_callback("GetShape", _GetShape)
	hc.add_callback("GetShape",ds,Rvir)
	hc.create(save_halos=True)
        #ad = ds.all_data()
    	#coordinatesDM = ad[("Halo","Coordinates")]
	#ellipsoid=a,b,c and rotation
	#add_callback(get_shape(coords,elipsoid,Rin,Rout))
	#ellipsoid=a,b,c and rotation
