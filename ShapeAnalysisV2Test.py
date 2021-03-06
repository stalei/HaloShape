#  © Shahram Talei @ 2019 The University of Alabama - All rights reserved.
#you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
# the main differnce in V2 is the COM instead of the center coordinates from rockstar

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
from yt.units import parsec, Msun, Mpc, cm
from mpl_toolkits.mplot3d import Axes3D

class MomentShape:
    def __init__(self,Rv,NBins):
        self.a=[0.]*NBins
        self.b=[0.]*NBins
        self.c=[0.]*NBins
        self.b_a=[0.]*NBins
        self.c_a=[0.]*NBins
        self.T=[0.]*NBins
        self.R=[0.]*NBins
    def SetBinShape(self,i,R,ellipsoid):
        ellipsoid.axis.sort()
        self.a[i]=ellipsoid.axis[0]
        self.b[i]=ellipsoid.axis[1]
        self.c[i]=ellipsoid.axis[2]
        self.b_a[i]=ellipsoid.b_a
        self.c_a[i]=ellipsoid.c_a
        self.T[i]=(ellipsoid.axis[0]**2.-ellipsoid.axis[1]**2.)/(ellipsoid.axis[0]**2.-ellipsoid.axis[2]**2.)
        self.R[i]=R
class EllipsoidShell:
    def __init__(self,axis,orientation,Rin,Rout):
        #I have to sort in this level
        aIndex=np.argmax(axis)
        cIndex=np.argmin(axis)
        if aIndex ==cIndex:
            self.axis=axis
            self.orientation=orientation
        else:
            bIndex=3-aIndex-cIndex # sum of indexes is always 3
            self.axis=np.array([axis[aIndex],axis[bIndex],axis[cIndex]])#.sort()
            self.orientation=(np.array([orientation[:,aIndex],orientation[:,bIndex],orientation[:,cIndex]])).T
        print(self.orientation)
        self.Rin=Rin
        self.Rout=Rout
        self.b_a=self.axis[1]/self.axis[0]
        self.c_a=self.axis[2]/self.axis[0]
    def IsInside(self,point,Rin,Rout):
        t=0
        pNew=np.array([0.0,0.0,0.0])
        for i in range(0,3):
            for j in range(0,3):
                pNew[i]+=point[j]*self.orientation[j,i] # each eig vec is [:,i] comp
        #for i in range(0,3):
        #    t+=(pNew[i]/(self.axis[i]))**2.
        #print("point, before-after projection")
        #print(point,pNew)
        t=np.sqrt(pNew[0]**2.+(pNew[1]/self.b_a)**2.+(pNew[2]/self.c_a)**2.)
        #print(t)
        if (t<=Rout and t>=Rin):
            return True
        else:
            return False
    def MomentTensor(self,coords,Rin,Rout):
        i=j=0
        c=0
        shape=[[0,0,0],[0,0,0],[0,0,0]]#np.array([[0,0,0,],[0,0,0,],[0,0,0]])
        s=0
        for i in range(0,3):
            for j in range(0,3):
                c=s=0
                print("s[%d,%d]"%(i,j))
                for point in coords:
                    #print("point before if:")
                    print(point)
                    if(self.IsInside(point,Rin,Rout)):
                        s+=point[i]*point[j]
                        c+=1
                #Mtot=c*m
                if c>0:
                    shape[i][j]=s/c
                print("particle count for %g<R<%g=%d"%(Rin,Rout,c))
        return shape
#def GetShellShape(self,sTen,)
def CompareShapes(ell1,ell2):
    AreSame=False;
    error=1.0e-1
    c=0
    #inP=np.array([0,0,0])
    dRatio=np.array([0,0])
    #we have to sort axis before we compare them
    for i in range(1,3):
        if(ell1.axis[0] !=0):
            dRatio[i-1]=np.abs((ell1.axis[i]/ell1.axis[0])/(ell2.axis[i]/ell2.axis[0]))
        #dRatio[i-1]=0
    #dori=np.abs(np.divide(ell1.orientations/ell2.orientations))
    d=dRatio[(dRatio<1+error) & (dRatio>1-error)]
    print(d)
    if(len(d)==2):
        for i in range(0,3):
            inP=np.array([0,0,0])
            for j in range(0,3):
                inP[i]+=ell1.orientation[j,i]*ell2.orientation[j,i]
            print(inP[i])
            if(inP[i]<((ell1.axis[i])**2.+error) and inP[i]>((ell1.axis[i])**2.-error)):
                c+=1
        if(c==3):
            AreSame=True
    return AreSame;



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
            #virialized
            Rcut=1.2*self.Rv
            pxh2=pxh[subr<Rcut]
            pyh2=pyh[subr<Rcut]
            pzh2=pzh[subr<Rcut]
            #with respect to center coordinate shift
            pxh2-=self.pos[0]
            pyh2-=self.pos[1]
            pzh2-=self.pos[2]
            coordsArray=np.array((pxh2,pyh2,pzh2))
            coordsHalo =coordsArray.T# we need to reshape
            return coordsHalo
        else:
            ad = snap.all_data()
            xp=ad[("all","particle_position_x")]
            yp=ad[("all","particle_position_y")]
            zp=ad[("all","particle_position_z")]
            print("this is x:")
            print(xp)
            coordsArray=np.array((xp,yp,zp))
            coordinatesDM =coordsArray.T# we need to reshape
            #coordinatesDM = ad[("Halo","Coordinates")]
            #IDsDM = ad[("Halo","ParticleIDs")]
            hid =h.id# halo.quantities['particle_identifier']
            #dds = halo.halo_catalog.data_ds
            #center = dds.arr([halo.quantities["particle_position_%s" % axis] \
            #for axis in "xyz"])
            c2=r2=0
            center =h.pos# dds.arr([halo.quantities["particle_position_%s" % axis] \
            Rcut=1.2*h.Rv
            coordsDM=np.array(coordinatesDM)#.v)
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
            print("r in extract:")
            print(r)
            #coords=coordinatesDM[np.abs(r-np.sqrt(c2)<Rvir)]
            coordsVir=coordsDM[r<Rcut]
            #IDsDmVir=IDsDM[r<Rcut]
            return coordsVir
    def ExtractShape(self,coords,NBins,NIteration):
        HShape=MomentShape(self.Rv,NBins)#define an empty shape object
        #let's create bins
        Rbins=np.linspace(0,self.Rv,NBins+1)
        Rs=[0]*NBins
        for i in range(0,NBins):
            Rs[i]=(Rbins[i]+Rbins[i+1])/2.
        print("bins:")
        print(Rbins)
        print("Rs:")
        print(Rs)
        #loop on bins starts here
        for i in range(0,NBins):
            convergence=False
            iteration=0
            #
            Rin=Rbins[i]
            Rout=Rbins[i+1]
            Rs[i]=(Rbins[i]+Rbins[i+1])/2.
            axis=np.array([Rout,Rout,Rout])
            orientation=np.identity(3)
            print(orientation)
            v0=reduce(mul,axis)
            #iteration starts here
            while(not(convergence)):
                print("Iteration:%d"%iteration)
                match=0
                ellshell=EllipsoidShell(axis,orientation,Rin,Rout)
                MTensor=ellshell.MomentTensor(coords,Rin,Rout)
                print(MTensor)
                axisNew,orientationNew=LA.eig(MTensor)
                v=reduce(mul,axisNew)
                norm=(v0/v)**(1./3.)
                axisNormalized=[ax*norm for ax in axisNew]
                print("axis, axisNormalized:")
                print(axisNew)
                print(axisNormalized)
                ellshellNew=EllipsoidShell(axisNormalized,orientationNew,Rin,Rout)
                if CompareShapes(ellshell,ellshellNew):
                    match+=1
                    print("we converged in shape for %f"%Rs[i])
                    convergence=True
                else:
                    match=0
                iteration+=1
                axis=axisNormalized
                orientation=orientationNew
                #if(match>=1):
                #    print("we converged in shape for %f"%Rs[i])
                #    convergence=True
                if (iteration>=NIteration):
                    convergence=True
                    print("The shape didn't converge but we reached the iteration limit for shell:%f"%Rs[i])
                if(convergence==True):
                    ellshape=EllipsoidShell(axis,orientation,Rin,Rout)
            HShape.SetBinShape(i,Rs[i],ellshape)
        return HShape


##################################################################
#tools for the shape



    #how to run: python ShapeAnalysis.py  bin_number iteration
    #example: $python ShapeAnalysisV2Test.py  5 30
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #parser.add_argument("snap", type=str)
    #parser.add_argument("halo",type=str)
    #parser.add_argument("PList",type=str)
    #parser.add_argument("UMass",type=float)
    #parser.add_argument("LMass",type=float)
    #parser.add_argument("contamination", type=int)
    #parser.add_argument("extractShape", type=int)
    parser.add_argument("NBins", type=int)
    parser.add_argument("NIteration", type=int)
    #parser.add_argument("excludeSubhalos", type=int)
    args = parser.parse_args()
    #read the snapshot file
    # in this test we create an arbitrary data instead of using an output file
    rRes=4
    tetRes=25
    fiRes=12
    n_particles = rRes*tetRes*fiRes
    Rvirial=0.2#200
    rpp=np.linspace(0.001,Rvirial-0.001,rRes)
    tetpp=np.linspace(-np.pi,np.pi,tetRes)
    fipp=np.linspace(0,2.*np.pi,fiRes)
    rMesh,tetMesh,fiMesh=np.meshgrid(rpp,tetpp,fipp)
    #rStrip=np.squeeze(np.array(rMesh))  x
    #rStrip=rMesh[:] x
    rStrip=np.reshape(rMesh,n_particles)
    tetStrip=np.reshape(tetMesh,n_particles)
    fiStrip=np.reshape(fiMesh,n_particles)
    #print(rMesh.shape)
    print(rStrip)
    ppx=rStrip*np.sin(tetStrip)*np.cos(fiStrip)
    ppy=rStrip*np.sin(tetStrip)*np.sin(fiStrip)
    ppz=rStrip*np.cos(tetStrip)
    fig0 = plt.figure(0,figsize=plt.figaspect(0.5))
    ax0 = fig0.add_subplot(121, projection='3d')
    ax0.scatter(ppx,ppy,ppz,c='black',alpha=0.9,marker='.',s=1)
    print("coords:")
    print(ppx)
    print(ppy)
    print(ppz)
    np.savetxt('x.out', ppx, delimiter=',')
    np.savetxt('y.out', ppy, delimiter=',')
    np.savetxt('z.out', ppz, delimiter=',')
    #ppx, ppy, ppz =1e2*np.random.normal(size=[3, n_particles])
    ppm = np.ones(n_particles)
    data = {'particle_position_x': ppx,'particle_position_y': ppy,'particle_position_z': ppz,'particle_mass': ppm}
    bbox = 1.1*np.array([[min(ppx), max(ppx)], [min(ppy), max(ppy)], [min(ppz), max(ppz)]])
    #snap = yt.load_particles(data, length_unit=Mpc, mass_unit=Msun, n_ref=256, bbox=bbox)
    snap = yt.load_particles(data, length_unit=cm, mass_unit=Msun, n_ref=256, bbox=bbox)
    #snap = yt.load(args.snap)#, unit_base=unit_base1)#,unit_system='galactic')
    #if snap is None:
    #    print ("Error, sorry, I couldn't read the snapshot.!")
    #    sys.exit(1)
    halos=0#np.genfromtxt(args.halo, skip_header=18)
    #if halos is None:
    #    print ("Error, sorry, I couldn't read the halo binary file.!")
    #    sys.exit(1)
    plist=0#np.genfromtxt(args.PList, skip_header=18,comments='#')
    #if plist is None:
    #    print ("Error, sorry, I couldn't read the halo binary file.!")
    #    sys.exit(1)
    #UpperMass=args.UMass
    #LowerMass=args.LMass
    #read the halo info
    #pnumh=np.array(halos[:,1])
    #MvH=np.array(halos[:,2])
    #RvH=np.array(halos[:,4])# in kpc
    #xH=np.array(halos[:,8])
    #yH=np.array(halos[:,9])
    #zH=np.array(halos[:,10])
    #IdH=np.array(halos[:,0])
    ph=n_particles#pnumh[(MvH>LowerMass) & (MvH<UpperMass)]
    Idh=0#IdH[(MvH>LowerMass) & (MvH<UpperMass)]
    Mvh=1.2e12#MvH[(MvH>LowerMass) & (MvH<UpperMass)]
    xh=0#xH[(MvH>LowerMass) & (MvH<UpperMass)]
    yh=0#yH[(MvH>LowerMass) & (MvH<UpperMass)]
    zh=0#zH[(MvH>LowerMass) & (MvH<UpperMass)]
    Rvh=Rvirial#RvH[(MvH>LowerMass) & (MvH<UpperMass)]
    #Rvh/=1000 # convert from kpc to Mpc
    halo=[]#*len(MvH)
    #print("found %d halos in the given interval"%len(Mvh))
    if True: #A #we already have a halo
    #    for i in range(0,len(Mvh)):
    #        halo.append(Halo(xh[i],yh[i],zh[i],Rvh[i],Mvh[i],ph[i],Idh[i]))
        halo.append(Halo(xh,yh,zh,Rvh,Mvh,ph,Idh))
        #if(args.contamination == 1):
        #    print("checking for contamination")
        #    for h in halo:
        #        h.contamination=h.IsContaminated(snap)
        print("ID ---- Mvirial ---- Rvirial ---- # of Ps ---- Contamination ---- [#p Vir ---- #p No Sub]")
        print("################################################################################")
        #if(args.extractShape==0):
        #    for h in halo:
        #        print("%d -- %g -- %f -- %d -- %d"%(h.id,h.Mv,h.Rv,h.pnum,h.contamination))
        #target=input("please select a halo id(-1 to quit):")
        #if target==-1:
        #    sys.exit(1)
        #print("let's extract the shape for:%d"%int(target))
        #tID=int(target)
        #print("####################################################")
        #start plots
        fig = plt.figure(1)
        fig.suptitle('Shapes')
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)
        ax1.legend(loc=2)
        ax1.title.set_text('b/a')
        ax2.title.set_text('c/a')
        ax3.title.set_text('T')
        ax4.title.set_text('c/a - b/a')
        if(True): # extract shape is already true and I don't want to change the structure
            #print("let's extract the shape")
            for h in halo:
                #print("%d -- %g -- %f -- %d -- %d -- %d -- %d"%(h.id,h.Mv,h.Rv,h.pnum,h.contamination,len(PcoordsSub),len(PcoordsNoSub)))
                #extract coordinates
                #
                PcoordsSub=h.ExtractParticles(snap,halos,plist,False)
                ax00 = fig0.add_subplot(122, projection='3d')
                ax00.scatter(PcoordsSub[:,0],PcoordsSub[:,1],PcoordsSub[:,2],c='black',alpha=0.9,marker='.',s=1)
                #PcoordsNoSub=h.ExtractParticles(snap,halos,plist,True)
                Shape=h.ExtractShape(PcoordsSub,args.NBins,args.NIteration)
                #ShapeNoSub=h.ExtractShape(PcoordsNoSub,args.NBins,args.NIteration)
                #print(Shape.R)
                #print(Shape.b_a)
                ax1.plot(Shape.R,Shape.b_a,linestyle='-',label=str(h.id))
                #ax1.plot(ShapeNoSub.R,ShapeNoSub.b_a,linestyle='-.',label="NoSub"+str(h.id))
                ax2.plot(Shape.R,Shape.c_a,linestyle='-',label=str(h.id))
                #ax2.plot(ShapeNoSub.R,ShapeNoSub.c_a,linestyle='-.',label="NoSub"+str(h.id))
                ax3.plot(Shape.R,Shape.T,linestyle='-',label=str(h.id))
                #ax3.plot(ShapeNoSub.R,ShapeNoSub.T,linestyle='-.',label="NoSub"+str(h.id))
                ax4.plot(Shape.c_a,Shape.b_a,'bo',label=str(h.id))
                #ax4.plot(ShapeNoSub.c_a,ShapeNoSub.b_a,'ro',label="NoSub"+str(h.id))
                print("%d -- %g -- %f -- %d -- %d -- %d -- %d"%(h.id,h.Mv,h.Rv,h.pnum,h.contamination,len(PcoordsSub),0))
        print("################################################################################")
    else: #A
        print("No halo found in the given range!")
        sys.exit(1)
    ax1.legend(loc=3)
    ax2.legend(loc=3)
    ax3.legend(loc=3)
    ax4.legend(loc=4)
    plt.show()
