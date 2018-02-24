
# read galaxia output
# generate input file for ubgaiaerror

import ebf
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from scipy import constants as const
from FortranFile import FortranFile
import struct

# output option 0: off, 1: on
# output UB Gaia error code ASCII data
flagubeasc=1
# output UB Gaia error code input binary data
flagube=1
# output ebf data
flagebf=0
# plot 
flagplot=1

# input filename
inputfile='../galaxia/galaxy1.ebf'
# ebf output filename
outputfile='galaxy1eq.ebf'

# magnitude limit
vmaglim=20.0

# reading the data
px=ebf.read(inputfile,'/px')
py=ebf.read(inputfile,'/py')
pz=ebf.read(inputfile,'/pz')
vy=ebf.read(inputfile,'/vy')
vx=ebf.read(inputfile,'/vx')
vz=ebf.read(inputfile,'/vz')
glon=ebf.read(inputfile,'/glon')
glat=ebf.read(inputfile,'/glat')
age=ebf.read(inputfile,'/Age')
feh=ebf.read(inputfile,'/FeH')
alfe=ebf.read(inputfile,'/Alpha')
smass=ebf.read(inputfile,'/Smass')
rad=ebf.read(inputfile,'/rad')
center=ebf.read(inputfile,'/Center')
teff=ebf.read(inputfile,'/Teff')
logg=ebf.read(inputfile,'/Grav')
lum=ebf.read(inputfile,'/lum')
mumag=ebf.read(inputfile,'/ubv_u')
mbmag=ebf.read(inputfile,'/ubv_b')
mvmag=ebf.read(inputfile,'/ubv_v')
mrmag=ebf.read(inputfile,'/ubv_r')
mimag=ebf.read(inputfile,'/ubv_i')
mjmag=ebf.read(inputfile,'/ubv_j')
mhmag=ebf.read(inputfile,'/ubv_h')
mkmag=ebf.read(inputfile,'/ubv_k')
exbvs=ebf.read(inputfile,'/exbv_schlegel')

# solar position and velocity
print 'Solar position assumed in Galaxia =',center[0],center[1],center[2]
print '      velocity =',center[3],center[4],center[5]

ns=len(mumag)

# apparent magnitude
vmag=mvmag+5.0*np.log10(rad*100.0)
imag=mimag+5.0*np.log10(rad*100.0)
# extinction
av=exbvs*3.315
ai=exbvs*1.940
# add extinction
vmage=vmag+av
image=imag+ai
vicole=vmage-image

# selecting apparent magnitude 
sindx=np.where(vmage < vmaglim)
ns=np.size(sindx)
print 'Ns(V<',vmaglim,')=',ns
pxs=px[sindx]
pys=py[sindx]
pzs=pz[sindx]
vxs=vx[sindx]
vys=vy[sindx]
vzs=vz[sindx]
glonas=glon[sindx]
glatas=glat[sindx]
vmages=vmage[sindx]
vicoles=vicole[sindx]

# adjusting the position
# no correction for the position
# pxs=pxs+center[0]
# pys=pys+center[1]
# pzs=pzs+center[2]
vxs=vxs+center[3]
vys=vys+center[4]
vzs=vzs+center[5]

# output ascii data for test
#f=open('testi-part.dat','w')
#i=0
#while i < ns:
# print >>f, "%f %f %f %f %f %f" %(pxs[i]-8.0,pys[i],pzs[i]
#  ,vxs[i],vys[i]+vrotmsun,vzs[i])
# i+=1
#f.close()

# 3 axes in Galactic coordinate
# North Celestial Pole (z-axis)
lzeq=np.radians(122.93193212)
bzeq=np.radians(27.12835496)
# RA,DEC=0,0 (x-axis)
lxeq=np.radians(96.33723825)
bxeq=np.radians(-60.18853909)
#  RA,DEC=90,0 (y-axis)
lyeq=np.radians(206.98916373)
byeq=np.radians(-11.42442440)
# transformation matrix (though it is array in python)
tmateqga=np.array([
  [np.cos(lxeq)*np.cos(bxeq),np.sin(lxeq)*np.cos(bxeq),np.sin(bxeq)],
  [np.cos(lyeq)*np.cos(byeq),np.sin(lyeq)*np.cos(byeq),np.sin(byeq)],
  [np.cos(lzeq)*np.cos(bzeq),np.sin(lzeq)*np.cos(bzeq),np.sin(bzeq)]
])

# transfer the coordinate: x-y plane: Galactic plane -> Cerestial equator

# position
# position in Galactic cartesian coordinate
posgas=np.vstack([pxs,pys,pzs])
# transfer to equatorial cartesian coordinate
poseqs=np.dot(tmateqga,posgas)

# velocity
# velocity in Galactic cartesian coordinate
velgas=np.vstack([vxs,vys,vzs])
# transfer to equatorial cartesian coordinate
veleqs=np.dot(tmateqga,velgas)


# R.A. Dec. 
radxys=np.sqrt(np.power(poseqs[0,:],2)+np.power(+poseqs[1,:],2))
rads=np.sqrt(np.power(poseqs[0,:],2)+np.power(poseqs[1,:],2)
  +np.power(poseqs[2,:],2))
# R.A.
alps=np.degrees(np.arccos(poseqs[0,:]/radxys))
i=0
while i<ns:
 if poseqs[1,i] < 0.0:
   alps[i]=360.0-alps[i]
 i+=1
dels=np.degrees(np.arcsin(poseqs[2,:]/rads))
equs=np.vstack([alps,dels])

# radial velocity
vrads=(veleqs[0,:]*poseqs[0,:]+veleqs[1,:]*poseqs[1,:]+veleqs[2,:]*poseqs[2,:])/rads
# proper motion
# mu_alpha (km/s)
valps=(-veleqs[0,:]*poseqs[1,:]+veleqs[1,:]*poseqs[0,:])/radxys
# mu_delta (km/s)
# vrad in x-y plane
vradxys=(veleqs[0,:]*poseqs[0,:]+veleqs[1,:]*poseqs[1,:])/radxys
vdels=(-vradxys*poseqs[2,:]+veleqs[2,:]*radxys)/rads

# constant for proper motion unit conversion
pmvconst=4.74047
# changing the unit from km/s to mas/yr
# km/s -> arcsec/y   vt=4.74 mu d(pc)
valps=1000.0*valps/pmvconst/(rads*1000.0)
vdels=1000.0*vdels/pmvconst/(rads*1000.0)

# test
# f=open('vposeq-test.dat','w')
# i=0
# while i < ns:
#  print >>f, "%f %f %f" %(px[i],py[i],pz[i])
#  print >>f, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f" %(poseqs[0,i],poseqs[1,i]
#   ,poseqs[2,i],veleqs[0,i],veleqs[1,i],veleqs[2,i]
#   ,alps[i],dels[i],vrads[i],valps[i]
#   ,vdels[i],rads[i],vxs[i],vys[i],vzs[i],pxs[i],pys[i],pzs[i])
#  i+=1
# f.close()


# stellar parameters for selected stars
ages=age[sindx]
# age Gyr
ages=np.power(10.0,age[sindx]-9.0)
fehs=feh[sindx]
alfes=alfe[sindx]
avs=av[sindx]
vicoles=vicole[sindx]
mvmags=mvmag[sindx]
smasss=smass[sindx]
teffs=teff[sindx]
# Teff K
teffs=np.power(10.0,teffs)
loggs=logg[sindx]

# bolometric magnitude
mbolsun=4.75
mbols=-2.5*lum+mbolsun

# distance in parsec
diss=rads*1000.0

# ASCII data for UB Fortran code
if flagubeasc:
  f=open('ubgaiae-in.asc','w')
  i=0
  while i < ns:
    print >>f, " %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e" %(alps[i],dels[i],diss[i],valps[i],vdels[i],vrads[i]
   ,teffs[i],logg[i],fehs[i],avs[i],vicoles[i],vmages[i],vmages[i]-vicoles[i]
   ,glatas[i],glonas[i])
    i+=1
  f.close()

# binary data for UB Fortran code
if flagube:
  f=FortranFile('ubgaiaein.bin',mode='w')
  nsarr=np.reshape(ns,1)
  f.writeInts(nsarr)
  i=0
  while i < ns:
    staro=np.array([alps[i],dels[i],diss[i],valps[i],vdels[i],vrads[i]
     ,teffs[i],logg[i],fehs[i],avs[i],vicoles[i],vmages[i]])
    f.writeReals(staro,prec='d')
    i+=1
  f.close()

# ebf output
if flagebf:
  ebf.write(outputfile,'/alpha',alps,'w')
  ebf.write(outputfile,'/delta',dels,'a')
  ebf.write(outputfile,'/distance',diss,'a')
  ebf.write(outputfile,'/mualpha',valps,'a')
  ebf.write(outputfile,'/mudelta',vdels,'a')
  ebf.write(outputfile,'/vrad',vdels,'a')
  ebf.write(outputfile,'/age',ages,'a')
  ebf.write(outputfile,'/feh',fehs,'a')
  ebf.write(outputfile,'/alfe',alfes,'a')
  ebf.write(outputfile,'/Av',avs,'a')
  ebf.write(outputfile,'/VIcole',vicoles,'a')
  ebf.write(outputfile,'/Mv',mvmags,'a')
  ebf.write(outputfile,'/Vape',vmages,'a')
  ebf.write(outputfile,'/smass',smasss,'a')
  ebf.write(outputfile,'/Teff',teffs,'a')
  ebf.write(outputfile,'/logg',loggs,'a')

if flagplot:
# test plot
# top panels
  gs1=gridspec.GridSpec(1,3)
  gs1.update(left=0.05,right=0.975,bottom=0.5,top=0.95,hspace=0.5)
# x-y plot
  plt.subplot(gs1[0],aspect='equal')
  plt.axis([-5.0,15.0,-10.0,10.0])
  plt.hexbin(pxs,pys,bins='log',gridsize=1000,cmap=cm.jet)
  plt.xlabel("x-y",fontsize=12,fontname="serif")
#plt.ylabel("y",fontsize=12,fontname="serif")
# x-z plot
  plt.subplot(gs1[1],aspect='equal')
  plt.axis([-5.0,15.0,-10.0,10.0])
  plt.hexbin(pxs,pzs,bins='log',gridsize=1000,cmap=cm.jet)
  plt.xlabel("x-z",fontsize=12,fontname="serif")
#plt.ylabel("z",fontsize=12,fontname="serif")
# Mv vs. V-I plot
  plt.subplot(gs1[2],aspect=0.2)
# hexbin plot
  plt.hexbin(vicoles,mvmags,bins='log',gridsize=200,cmap=cm.jet)
  plt.axis([-0.5,3.0,15.0,-10.0])
# labes
  plt.xlabel(r"$\rm M_V vs. V-I$",fontsize=12,fontname="serif")
#plt.ylabel(r"$\rm M_V$",fontsize=12,fontname="serif")

# bottom panel
  gs2=gridspec.GridSpec(1,2)
  gs2.update(left=0.1,right=0.975,bottom=0.05,top=0.4)
# Galactic coordinate
  plt.subplot(gs2[0],aspect='equal')
  plt.hexbin(glonas,glatas,bins='log',gridsize=100,cmap=cm.jet)
  plt.axis([0.0,360.0,-90.0,90.0])
  plt.xlabel("l",fontsize=12,fontname="serif")
  plt.ylabel("b",fontsize=12,fontname="serif")
# equatorial coordinate
  plt.subplot(gs2[1],aspect='equal')
  plt.hexbin(equs[0,:],equs[1,:],bins='log',gridsize=100,cmap=cm.jet)
  plt.axis([0.0,360.0,-90.0,90.0])
  plt.xlabel(r"$\alpha$",fontsize=12,fontname="serif")
  plt.ylabel(r"$\delta$",fontsize=12,fontname="serif")
#cb=plt.colorbar()
  plt.show()
  plt.savefig('galaxia-stars.png')
