#  © Shahram Talei @ 2019 The University of Alabama - All rights reserved.
#you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.legend_handler import HandlerLine2D
from mpl_toolkits.mplot3d import Axes3D
import argparse
plt.rcParams["font.size"] =12

def PlotShpere(ax,R,xc,yc,zc):
	u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
	x = xc+(R/1000.)*np.cos(u)*np.sin(v)
	y = yc+(R/1000.)*np.sin(u)*np.sin(v)
	z = zc+(R/1000.)*np.cos(v)
	ax.plot_wireframe(x, y, z, color="gray",alpha=0.3)
########################################################### Data source
# #id num_p mvir mbound_vir rvir vmax rvmax vrms x y z vx vy vz Jx Jy Jz E Spin PosUncertainty VelUncertainty
#bulk_vx bulk_vy bulk_vz BulkVelUnc n_core m200b m200c m500c m2500c Xoff Voff spin_bullock b_to_a c_to_a A[x] A[y] A[z]
#b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c) A[z](500c) Rs Rs_Klypin T/|U| M_pe_Behroozi M_pe_Diemer
#Halfmass_Radius idx i_so i_ph num_cp mmetric

#how to run: python SelectHalo.py halo_catalog num_limit M_high M_low Plot_dwarfs
#example: $python SelectHalo.py halos_0.0_G.ascii 250 1.25e12 1.2e12 1

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("HaloCatalog",type=str)
	parser.add_argument("num_limit", type=int)
	parser.add_argument("M_high", type=float)
	parser.add_argument("M_low", type=float)
	parser.add_argument("plotDwarfs", type=int)
	args = parser.parse_args()
	data=np.genfromtxt(args.HaloCatalog, skip_header=18)#,names=True, skip_header=5)
	IdAll=np.array(data[:,0])
	#CountAll= len(id)
	NumPAll=np.array(data[:,1])
	MvirAll=np.array(data[:,2])
	RvirAll=np.array(data[:,4])
	XAll=np.array(data[:,8])
	YAll=np.array(data[:,9])
	ZAll=np.array(data[:,10])
	## Criterias

	MaxMass=args.M_high #1.2e12
	MinMass=args.M_low#0.98e11
	NumLimit=args.num_limit#250


	#to see smaller satellites
	IdDwarfs=IdAll[MvirAll<MinMass]
	RvirDwarfs=RvirAll[MvirAll<MinMass]
	XDwarfs=XAll[MvirAll<(MinMass)]
	YDwarfs=YAll[MvirAll<(MinMass)]
	ZDwarfs=ZAll[MvirAll<(MinMass)]


	# remove hlos with particle number less than NumLimit
	IdHigh=IdAll[NumPAll>NumLimit]
	MvirHigh=MvirAll[NumPAll>NumLimit]
	RvirHigh=RvirAll[NumPAll>NumLimit]
	XHigh=XAll[NumPAll>NumLimit]
	YHigh=YAll[NumPAll>NumLimit]
	ZHigh=ZAll[NumPAll>NumLimit]


	# remove masses smaller than MinMass
	IdAbove=IdHigh[MvirHigh>MinMass]
	MvirAbove=MvirHigh[MvirHigh>MinMass]
	RvirAbove=RvirHigh[MvirHigh>MinMass]
	XAbove=XHigh[MvirHigh>MinMass]
	YAbove=YHigh[MvirHigh>MinMass]
	ZAbove=ZHigh[MvirHigh>MinMass]


	# remove masses greater than MaxMass
	Id=IdAbove[MvirAbove<MaxMass]
	Mvir=MvirAbove[MvirAbove<MaxMass]
	Rvir=RvirAbove[MvirAbove<MaxMass]
	X=XAbove[MvirAbove<MaxMass]
	Y=YAbove[MvirAbove<MaxMass]
	Z=ZAbove[MvirAbove<MaxMass]



	#print(XAll[IdAll==33678],YAll[IdAll==33678],ZAll[IdAll==33678])
	#print(XAll[IdAll==33743],YAll[IdAll==33743],ZAll[IdAll==33743])


	print("Min # of particles in a halo is:", NumLimit)

	print("No of halos within mass limit: ", "%.4g"%MinMass, "", "%.4g"%MaxMass,"is:",len(Id))



	########################################################## plots


	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	label =str(Id)


	for i in range(len(X)):
		PlotShpere(ax,Rvir[i],X[i],Y[i],Z[i])
		ax.scatter(X[i],Y[i],Z[i],c='black', alpha=0.9, marker='.',s=15)#Rvir[i])#15)#,s= (10.0*np.log10(Mvir[i])))
		ax.text(X[i],Y[i],Z[i],'%s' %str(Id[i]), size=7)
	if args.plotDwarfs==1:
		ax.scatter(XDwarfs,YDwarfs,ZDwarfs,c='r',alpha=0.8,marker='o',s=8)

	fig.set_size_inches(14,8)
	ax.set_xlabel('X (Mpc)')
	ax.set_ylabel('Y (Mpc)')
	ax.set_zlabel('Z (Mpc)')

	plt.show()
