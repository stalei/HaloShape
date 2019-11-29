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

class Halo:
    def __init__(self,x,y,z,R,M,pnum,id):
        center=np.zeros(3)
        center[0]=x
        center[1]=y
        center[2]=z
        self.pos=center
        self.Rv=R
        self.Mv=M
        self.pnum=pnum
        self.id=id
        self.contamination=0
    #look for low resolution particle contamination
    def IsContaminated(self,ds):#snap):
        center =self.pos# dds.arr([halo.quantities["particle_position_%s" % axis] \
        Rv=self.Rv
        count=0
        boundary=np.zeros((2,3))
        ad = ds.all_data()
        coordinatesDM = ad[("Halo","Coordinates")]
        coordinatesLowRes = ad[("Bndry","Coordinates")]
        DMCount=len(coordinatesDM)
        LowResCount=len(coordinatesLowRes)
        print("number of High resolution particles:%d"%DMCount)
        print("number of Low resolution particles:%d"%LowResCount)
        for i in range(0,3):
            boundary[0,i]=np.min(coordinatesDM[:,i]) # (min max, x y z)
            boundary[1,i]=np.max(coordinatesDM[:,i])
        print("high resolution particles are in box:")
        print(boundary)
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
        r2=((xLR4.v-center[0])**2.+(yLR4.v-center[1])**2.+(zLR4.v-center[2])**2.)
        r=np.sqrt(r2)
        rContamination=r[r<(Rv)]
        count=len(rContamination)
        print("Halo %d at:" % (h.id))
        print(center)
        print("is contaminated with: %d low resolution particles at δr="%count)
        print(rContamination)
        return count
    def ExtractParticles(self,snap,halos,plist,ExcludeSubs):
        if ExcludeSubs==True:
            #
            subPs=plist#np.genfromtxt(plist, skip_header=18,comments='#')
            px=np.array(subPs[:,0])
            py=np.array(subPs[:,1])
            pz=np.array(subPs[:,2])
            pHids=np.array(subPs[:,9])
            pxh=px[pHids==self.id]
            pyh=py[pHids==self.id]
            pzh=pz[pHids==self.id]
            subr=np.sqrt((pxh-self.pos[0])**2.+(pyh-self.pos[1])**2.+(pzh-self.pos[2])**2.)
            pxh2=pxh[subr<self.Rv]
            pyh2=pyh[subr<self.Rv]
            pzh2=pzh[subr<self.Rv]
            coordsArray=np.array((pxh2,pyh2,pzh2))
            coordsHalo =coordsArray.T# we need to reshape
            return coordsHalo
        else:
            ad = snap.all_data()
            coordinatesDM = ad[("Halo","Coordinates")]
            IDsDM = ad[("Halo","ParticleIDs")]
            hid =h.id# halo.quantities['particle_identifier']
            #dds = halo.halo_catalog.data_ds
            #center = dds.arr([halo.quantities["particle_position_%s" % axis] \
            #for axis in "xyz"])
            c2=r2=0
            #center*=3.24078e-25
            #center=center.in_units('Mpc')
            center =h.pos# dds.arr([halo.quantities["particle_position_%s" % axis] \
            Rvir=h.Rv
            #print("halo center (Mpc):")
            #print(center)
            coordsDM=np.array(coordinatesDM.v)
            for i in range(0,3):
                #print(coordsDM[:,i])
                #print(center[i])
                coordsDM[:,i]-=center[i]
                #print(coordsDM[:,i])
                #print(i)
            #print(coordinatesDM)
            for i in range(0,3):
                r2+=(coordsDM[:,i])**2.
                #c2+=(center[i])**2.
            r=np.sqrt(r2)
            #print(r)
            #coords=coordinatesDM[np.abs(r-np.sqrt(c2)<Rvir)]
            coordsVir=coordsDM[r<Rvir]
            IDsDmVir=IDsDM[r<Rvir]
            return coordsVir

##################################################################
#tools for the shape



    #how to run: python ShapeAnalysis.py snapshot_file halo_catalog particles_list check_contamination extract_shape bin_number iteration
    #example: $python ShapeAnalysisV1.py snap_264 halos_0.0.ascii halos_0.0.particles 1.4e12 1.1e12   1 1 5 3
    #if don't want to exclude subhalo particles, simply pass 0 for particle list file
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("snap", type=str)
    parser.add_argument("halo",type=str)
    parser.add_argument("PList",type=str)
    parser.add_argument("UMass",type=float)
    parser.add_argument("LMass",type=float)
    parser.add_argument("contamination", type=int)
    parser.add_argument("extractShape", type=int)
    parser.add_argument("BinNum", type=int)
    parser.add_argument("IterationLim", type=int)
    #parser.add_argument("excludeSubhalos", type=int)
    args = parser.parse_args()
    #read the snapshot file
    snap = yt.load(args.snap)#, unit_base=unit_base1)#,unit_system='galactic')
    if snap is None:
        print ("Error, sorry, I couldn't read the snapshot.!")
        sys.exit(1)
    halos=np.genfromtxt(args.halo, skip_header=18)
    if halos is None:
        print ("Error, sorry, I couldn't read the halo binary file.!")
        sys.exit(1)
    plist=np.genfromtxt(args.PList, skip_header=18,comments='#')
    if plist is None:
        print ("Error, sorry, I couldn't read the halo binary file.!")
        sys.exit(1)
    UpperMass=args.UMass
    LowerMass=args.LMass
    #read the halo info
    pnumh=np.array(halos[:,1])
    MvH=np.array(halos[:,2])
    RvH=np.array(halos[:,4])# in kpc
    xH=np.array(halos[:,8])
    yH=np.array(halos[:,9])
    zH=np.array(halos[:,10])
    IdH=np.array(halos[:,0])
    ph=pnumh[(MvH>LowerMass) & (MvH<UpperMass)]
    Idh=IdH[(MvH>LowerMass) & (MvH<UpperMass)]
    Mvh=MvH[(MvH>LowerMass) & (MvH<UpperMass)]
    xh=xH[(MvH>LowerMass) & (MvH<UpperMass)]
    yh=yH[(MvH>LowerMass) & (MvH<UpperMass)]
    zh=zH[(MvH>LowerMass) & (MvH<UpperMass)]
    Rvh=RvH[(MvH>LowerMass) & (MvH<UpperMass)]
    Rvh/=1000 # convert from kpc to Mpc
    halo=[]#*len(MvH)
    print("found %d halos in the given interval"%len(Mvh))
    if len(Mvh)>0:

        for i in range(0,len(Mvh)):
            halo.append(Halo(xh[i],yh[i],zh[i],Rvh[i],Mvh[i],ph[i],Idh[i]))
        if(args.contamination == 1):
            print("checking for contamination")
            for h in halo:
                h.contamination=h.IsContaminated(snap)
        print("ID ---- Mvirial ---- Rvirial ---- # of Ps ---- Contamination ---- [#p Vir ---- #p No Sub]")
        print("################################################################################")
        if(args.extractShape==0):
            for h in halo:
                print("%d -- %g -- %f -- %d -- %d"%(h.id,h.Mv,h.Rv,h.pnum,h.contamination))
        #target=input("please select a halo id(-1 to quit):")
        #if target==-1:
        #    sys.exit(1)
        #print("let's extract the shape for:%d"%int(target))
        #tID=int(target)
        #print("####################################################")
        if(args.extractShape==1):
            #print("let's extract the shape")
            for h in halo:
                #extract coordinates
                #
                PcoordsSub=h.ExtractParticles(snap,halos,plist,False)
                PcoordsNoSub=h.ExtractParticles(snap,halos,plist,True)
                print("%d -- %g -- %f -- %d -- %d -- %d -- %d"%(h.id,h.Mv,h.Rv,h.pnum,h.contamination,len(PcoordsSub),len(PcoordsNoSub)))
        print("################################################################################")
