
# To plot distance error vs. distance   
# Errors are calcualted by a simple error propagation.

import ebf
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import constants as const
from FortranFile import FortranFile
import struct

# input parameters
# for distance vs. error
# number of bin
nbin=12
# number of rmax (kpc)
rmax=100.0

# input file name
inputfileint='gaiaei-out.bin'
inputfiledb='gaiaed-out.bin'

# reading ASCIIthe data
# rdata=np.loadtxt('ubgaiae-out.dat')
# reading binary data
f=FortranFile(inputfileint)
nset=f.readInts()
print ' input Number of stars =',nset
f.close()
f=FortranFile(inputfiledb)
i=0
rdata=np.zeros((nset,45))
while i < nset:
  rdata[i,:]=f.readReals('d')
  i+=1
f.close()

print ' Number of stars=',nset
# no error data
alpha=rdata[:,0]
delta=rdata[:,1]
parad=rdata[:,2]
mua=rdata[:,3]
mud=rdata[:,4]
vrad=rdata[:,5]
# error
alphaerr=rdata[:,12]
deltaerr=rdata[:,13]
paraderr=rdata[:,14]
muaerr=rdata[:,15]
muderr=rdata[:,16]
vraderr=rdata[:,17]

# only chose positive parallax and < 100 kpc away
#sindx=np.where(paradet > 1.0e-5)
# + parallax error < 100 mus
#sindx=np.where(np.logical_and(paradet > 1.0e-4,
# paradeterr/paradet < 0.1))

# get Cartician position in equatorial coordidate
# mas -> distance in kpc
dist=1.0/parad
# x,y,z position from the sun's position
rxy=dist*np.cos(delta)
px=rxy*np.cos(alpha)
py=rxy*np.sin(alpha)
pz=dist*np.sin(delta)
# vx,vy,vz
# arcsec/yr -> km/s
# constant for proper motion unit conversion
pmvconst=4.74047
valpha=pmvconst*mua/parad
vdelta=pmvconst*mud/parad
vxy=vrad*np.cos(delta)-vdelta*np.sin(delta)
vx=vxy*np.cos(alpha)-valpha*np.sin(alpha)
vy=vxy*np.sin(alpha)+valpha*np.cos(alpha)
vz=vrad*np.sin(delta)+vdelta*np.cos(delta)
# error
# distance mas -> kpc
#disterr=1.0/(parad+paraderr)
#disterr=dist-disterr
# error propagation
#disterr=np.log10(np.sqrt((1.0/np.power(parad,4))*np.power(paraderr,2)))
disterr=np.sqrt((1.0/np.power(parad,4))*np.power(paraderr,2))
# velocity
valphaerr=pmvconst*muaerr/parad
vdeltaerr=pmvconst*muderr/parad

# output ascii data for test
# f=open('detest.dat','w')
# i=0
# while i < nset:
# print >>f, "%f %f %f %f %f %f %f %f %f %f" %(px[i],py[i],pz[i]
#  ,vx[i],vy[i],vz[i],dist[i],disterr[i]
#  ,parad[i],paraderr[i])
# i+=1
# f.close()

# select 
sindx=np.where(np.logical_and(dist < 15.0
 ,np.logical_and(disterr < 15.0,disterr > 0.001)))
# selected particles for plot
#
pdist=dist[sindx]
pdisterr=disterr[sindx]
#
print ' nsplot=',len(pdist)

plt.subplot(111)
# labes
plt.xlabel("Distance (kpc)",fontsize=18,fontname="serif")
plt.ylabel("Distance Error (kpc)",fontsize=18,fontname="serif",style="normal")
# hexbin plot
plt.hexbin(pdist,pdisterr,yscale='log',bins='log',gridsize=80,cmap=cm.jet)
# plot mean
#plt.errorbar(rmean_r,vrotm_r,yerr=vrotsig_r,fmt='ow')
#plt.axis([pxmin,pxmax,pymin,pymax])
#plt.axis([0.0,15.0,-3.0,1.2])
plt.xticks(fontsize=18)
plt.yticks([0.01,0.1,1.0,10.0],("0.01","0.1","1.0","10.0"),fontsize=18)
cb=plt.colorbar()
cb.set_label('Log (Counts/100)')

#
plt.show()
plt.savefig('disterr.png')

