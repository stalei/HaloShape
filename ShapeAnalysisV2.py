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

import csv

class MomentShape:
    def __init__(self,Rv,NBins):
        self.a=[0.]*NBins
        self.b=[0.]*NBins
        self.c=[0.]*NBins
        self.A=[0.,0.,0.]*NBins
        self.B=[0.,0.,0.]*NBins
        self.C=[0.,0.,0.]*NBins
        self.Aunit=[0.,0.,0.]*NBins
        self.Bunit=[0.,0.,0.]*NBins
        self.Cunit=[0.,0.,0.]*NBins
        self.b_a=[0.]*NBins
        self.c_a=[0.]*NBins
        self.T=[0.]*NBins
        self.R=[0.]*NBins
    def SetBinShape(self,i,R,ellipsoid):
        ellipsoid.axis.sort()
        self.a[i]=ellipsoid.axis[0]
        self.b[i]=ellipsoid.axis[1]
        self.c[i]=ellipsoid.axis[2]
        self.A[i]=ellipsoid.orientation[:,0]
        self.B[i]=ellipsoid.orientation[:,1]
        self.C[i]=ellipsoid.orientation[:,2]
        self.Aunit[i] = self.A[i] / (self.A[i]**2).sum()**0.5
        self.Bunit[i] = self.B[i] / (self.B[i]**2).sum()**0.5
        self.Cunit[i] = self.C[i] / (self.C[i]**2).sum()**0.5
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
        t=np.sqrt(pNew[0]**2.+(pNew[1]/self.b_a)**2.+(pNew[2]/self.c_a)**2.)
        #print(t)
        if (t<=Rout and t>=Rin):
            return True
        else:
            return False
    def MomentTensor(self,coords,Rin,Rout):
        i=j=0
        c=0
        shape=[[1,0,0],[0,1,0],[0,0,1]]#np.array([[0,0,0,],[0,0,0,],[0,0,0]])
        s=0
        for i in range(0,3):
            for j in range(0,3): #it is symmetric and we can do range(i,3) to get it faster and later s.T+s and s[ii]/=2
                c=s=0
                print("s[%d,%d]"%(i,j))
                for point in coords:
                    #print("point before if:")
                    #print(point)
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
    def __init__(self,x,y,z,R,M,pnum,id,jx,jy,jz):
        center=np.zeros(3)
        am=np.zeros(3)
        center[0]=x
        center[1]=y
        center[2]=z
        am[0]=jx
        am[1]=jy
        am[2]=jz
        self.pos=center
        self.Rv=R
        self.Mv=M
        self.pnum=pnum
        self.id=id
        self.contamination=0
        self.AngMom=am
        self.AngMomUnit=am/np.sqrt(jx**2.+jy**2.+jz**2.)
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
            coordinatesDM = ad[("Halo","Coordinates")]
            IDsDM = ad[("Halo","ParticleIDs")]
            hid =h.id# halo.quantities['particle_identifier']
            #dds = halo.halo_catalog.data_ds
            #center = dds.arr([halo.quantities["particle_position_%s" % axis] \
            #for axis in "xyz"])
            c2=r2=0
            center =h.pos# dds.arr([halo.quantities["particle_position_%s" % axis] \
            Rcut=1.2*h.Rv
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
            coordsVir=coordsDM[r<Rcut]
            IDsDmVir=IDsDM[r<Rcut]
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



    #how to run: python ShapeAnalysis.py snapshot_file halo_catalog particles_list check_contamination extract_shape bin_number iteration
    #example: $python ShapeAnalysisV1.py snap_264 halos_0.0.ascii halos_0.0.particles 1.4e12 1.1e12   1 1 5 3
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("snap", type=str)
    parser.add_argument("halo",type=str)
    parser.add_argument("PList",type=str)
    parser.add_argument("UMass",type=float)
    parser.add_argument("LMass",type=float)
    parser.add_argument("contamination", type=int)
    parser.add_argument("extractShape", type=int)
    parser.add_argument("NBins", type=int)
    parser.add_argument("NIteration", type=int)
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
    JxH=np.array(halos[:,14])
    JyH=np.array(halos[:,15])
    JzH=np.array(halos[:,16])
    ph=pnumh[(MvH>LowerMass) & (MvH<UpperMass)]
    Idh=IdH[(MvH>LowerMass) & (MvH<UpperMass)]
    Mvh=MvH[(MvH>LowerMass) & (MvH<UpperMass)]
    xh=xH[(MvH>LowerMass) & (MvH<UpperMass)]
    yh=yH[(MvH>LowerMass) & (MvH<UpperMass)]
    zh=zH[(MvH>LowerMass) & (MvH<UpperMass)]
    Rvh=RvH[(MvH>LowerMass) & (MvH<UpperMass)]
    Jxh=JxH[(MvH>LowerMass) & (MvH<UpperMass)]
    Jyh=JyH[(MvH>LowerMass) & (MvH<UpperMass)]
    Jzh=JzH[(MvH>LowerMass) & (MvH<UpperMass)]
    Rvh/=1000 # convert from kpc to Mpc
    halo=[]#*len(MvH)
    print("found %d halos in the given interval"%len(Mvh))
    if len(Mvh)>0: #A

        for i in range(0,len(Mvh)):
            halo.append(Halo(xh[i],yh[i],zh[i],Rvh[i],Mvh[i],ph[i],Idh[i],Jxh[i],Jyh[i],Jzh[i]))
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
        #start plots
        fig = plt.figure(1)
        fig.suptitle('Shapes')
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)
        ax1.legend(loc=2)
        ax1.title.set_text('b/a')
        ax1.set_xlabel('R (Mpc)')
        ax1.set_ylabel('b/a')
        ax2.title.set_text('c/a')
        ax2.set_xlabel('R (Mpc)')
        ax2.set_ylabel('c/a')
        ax3.title.set_text('T')
        ax3.set_xlabel('R (Mpc)')
        ax3.set_ylabel('T')
        ax4.title.set_text('c/a - b/a')
        ax4.set_xlabel('b/a')
        ax4.set_ylabel('c/a')
        #2nd fig
        fig2=plt.figure(2)
        fig2.suptitle('Orientations')
        ax21=fig2.add_subplot(121)
        ax21.title.set_text('residual')
        ax21.set_xlabel('R (Mpc)')
        ax21.set_ylabel('inner product')
        ax22=fig2.add_subplot(122)
        ax22.title.set_text('halo angular momentum')
        ax22.set_xlabel('R (Mpc)')
        ax22.set_ylabel('inner product')
        if(args.extractShape==1):
            #print("let's extract the shape")
            for h in halo:
                #print("%d -- %g -- %f -- %d -- %d -- %d -- %d"%(h.id,h.Mv,h.Rv,h.pnum,h.contamination,len(PcoordsSub),len(PcoordsNoSub)))
                #extract coordinates
                #
                PcoordsSub=h.ExtractParticles(snap,halos,plist,False)
                PcoordsNoSub=h.ExtractParticles(snap,halos,plist,True)
                Shape=h.ExtractShape(PcoordsSub,args.NBins,args.NIteration)
                ShapeNoSub=h.ExtractShape(PcoordsNoSub,args.NBins,args.NIteration)
                #print(Shape.R)
                #print(Shape.b_a)
                ax1.plot(Shape.R,Shape.b_a,linestyle='-',label=str(h.id))
                ax1.plot(ShapeNoSub.R,ShapeNoSub.b_a,linestyle='-.',label="NoSub"+str(h.id))
                ax2.plot(Shape.R,Shape.c_a,linestyle='-',label=str(h.id))
                ax2.plot(ShapeNoSub.R,ShapeNoSub.c_a,linestyle='-.',label="NoSub"+str(h.id))
                ax3.plot(Shape.R,Shape.T,linestyle='-',label=str(h.id))
                ax3.plot(ShapeNoSub.R,ShapeNoSub.T,linestyle='-.',label="NoSub"+str(h.id))
                ax4.plot(Shape.c_a,Shape.b_a,'bo',label=str(h.id))
                ax4.plot(ShapeNoSub.c_a,ShapeNoSub.b_a,'ro',label="NoSub"+str(h.id))
                #orientation plots, I can add these in SetBinShape
                InnProA=[0]*args.NBins
                InnProResA=[0]*args.NBins
                InnProB=[0]*args.NBins
                InnProResB=[0]*args.NBins
                InnProC=[0]*args.NBins
                InnProResC=[0]*args.NBins
                print(Shape.Aunit)
                aa0=np.array(Shape.Aunit[0])
                bb0=np.array(Shape.Bunit[0])
                cc0=np.array(Shape.Cunit[0])
                print('this is aa0:')
                print(aa0)
                #index=0
                for index in range(0,args.NBins):
                    suma=0
                    sumResa=0
                    sumb=0
                    sumResb=0
                    sumc=0
                    sumResc=0
                    aa=np.array(Shape.Aunit[index])
                    bb=np.array(Shape.Bunit[index])
                    cc=np.array(Shape.Cunit[index])
                    for k in range(0,3):
                        suma+=aa[k]*h.AngMomUnit[k]
                        sumResa+=aa[k]*aa0[k]
                        sumb+=bb[k]*h.AngMomUnit[k]
                        sumResb+=bb[k]*bb0[k]
                        sumc+=cc[k]*h.AngMomUnit[k]
                        sumResc+=cc[k]*cc0[k]
                    InnProA[index]=suma
                    InnProResA[index]=sumResa
                    InnProB[index]=sumb
                    InnProResB[index]=sumResb
                    InnProC[index]=sumc
                    InnProResC[index]=sumResc
                    #index+=1
                ax21.plot(Shape.R,InnProA,linestyle='-',label="A-"+str(h.id))
                ax22.plot(Shape.R,InnProResA,linestyle='-',label="A-"+str(h.id))
                ax21.plot(Shape.R,InnProB,linestyle='-.',label="B-"+str(h.id))
                ax22.plot(Shape.R,InnProResB,linestyle='-.',label="B-"+str(h.id))
                ax21.plot(Shape.R,InnProC,linestyle=':',label="C-"+str(h.id))
                ax22.plot(Shape.R,InnProResC,linestyle=':',label="C-"+str(h.id))
                # now let's save the file
                #first we need a header with the halo info
                #
                header=str(h.id)+str(h.AngMom)+str(h.Mv)+str(h.Rv)
                #2nd we need to put all useful info in a matrix and write in a csv file
                OutputData=(np.array([Shape.R,Shape.b_a,Shape.c_a,Shape.T,Shape.A,Shape.B,Shape.C,InnProA,InnProResA,InnProB,InnProResB,InnProC,InnProResC])).T
                fName=str(h.id)+".csv"
                with open(fName, 'w') as f:
                    w = csv.writer(f, delimiter=';')
                    w.writerow(header)
                    #for row in OutputData:#zip(l1, l2, (str(x)+'%' for x in l3)):
                    #    w.writerow(row)
                    for row in range(0,args.NBins):
                        Record=np.array([Shape.R[row],Shape.b_a[row],Shape.c_a[row],Shape.T[row],InnProA[row],InnProResA[row],InnProB[row],InnProResB[row],InnProC[row],InnProResC[row]])
                        w.writerow(Record)
                #print(Shape.A)
                #print(InnPro)
                #print(InnProRes)
                print("%d -- %g -- %f -- %d -- %d -- %d -- %d"%(h.id,h.Mv,h.Rv,h.pnum,h.contamination,len(PcoordsSub),len(PcoordsNoSub)))
        print("################################################################################")
    else: #A
        print("No halo found in the given range!")
        sys.exit(1)
    ax1.legend(loc=3)
    ax2.legend(loc=3)
    ax3.legend(loc=3)
    ax4.legend(loc=4)
    ax21.legend(loc=2)
    ax22.legend(loc=4)
    plt.show()
