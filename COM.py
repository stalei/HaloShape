#  © Shahram Talei @ 2021 The University of Alabama - All rights reserved.
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

plt.rcParams["font.size"] =12

    #how to run: python COM.py snapshot_file
    #example: $python COM.py snap_264
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("snap", type=str)
    args = parser.parse_args()
    snap = yt.load(args.snap)
    ad = snap.all_data()
    coordinatesDM = ad[("Halo","Coordinates")]
    p=np.array(coordinatesDM)
    #print(p[:,1].shape)
    #print(p[1:])
    x0=50.9522
    y0=53.1411
    z0=47.4905
    Rv=0.15 #200kpc
    px=p[:,0]
    py=p[:,1]
    pz=p[:,2]
    dx=px-x0
    dy=py-y0
    dz=pz-z0
    r2=dx*dx+dy*dy+dz*dz
    r=np.sqrt(r2)
    x=px[r<Rv]
    y=py[r<Rv]
    z=pz[r<Rv]
    print(len(x))
    com_x=np.sum(x)/len(x)
    com_y=np.sum(y)/len(y)
    com_z=np.sum(z)/len(z)
    print("COM (x,y,z):%g,%g,%g"%(com_x,com_y,com_z))
