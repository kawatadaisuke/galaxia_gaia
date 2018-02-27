
# convert gaia error data to fits file

import pyfits
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import constants as const
from FortranFile import FortranFile
import struct
from galpy.util import bovy_coords

# input parameters
# for distance vs. error
# number of bin
nbin=12
# number of rmax (kpc)
rmax=100.0
# magnitude limit
vmaglim=20.0

# input file name
inputfileint='gaiaei-out.bin'
inputfiledb='gaiaed-out.bin'
# output file
outfile='galaxia_gaia.fits'

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
# radian -> degree
rdata[:,0]=180.0*rdata[:,0]/np.pi
rdata[:,1]=180.0*rdata[:,1]/np.pi
rdata[:,6]=180.0*rdata[:,6]/np.pi
rdata[:,7]=180.0*rdata[:,7]/np.pi
rdata[:,12]=180.0*rdata[:,12]/np.pi
rdata[:,13]=180.0*rdata[:,13]/np.pi

# Galactic coordinates
# True l and b
RA_true=rdata[:,0]
DEC_true=rdata[:,1]
Tllbb=bovy_coords.radec_to_lb(RA_true,DEC_true,degree=True,epoch=2000.0)
GLON_true=Tllbb[:,0]
GLAT_true=Tllbb[:,1]
# True pmGLON, pmGLAT
pmRA_true=rdata[:,3]
pmDEC_true=rdata[:,4]
Tpmllbb=bovy_coords.pmrapmdec_to_pmllpmbb(pmRA_true,pmDEC_true \
   ,RA_true,DEC_true,degree=True,epoch=2000.0)
pmGLON_true=Tpmllbb[:,0]
pmGLAT_true=Tpmllbb[:,1]
# observed 
RA_obs=rdata[:,6]
DEC_obs=rdata[:,7]
pmRA_obs=rdata[:,9]
pmDEC_obs=rdata[:,10]
Tpmllbb=bovy_coords.pmrapmdec_to_pmllpmbb(pmRA_obs,pmDEC_obs \
   ,RA_true,DEC_true,degree=True,epoch=2000.0)
pmGLON_obs=Tpmllbb[:,0]
pmGLAT_obs=Tpmllbb[:,1]
# error
e_pmRA=rdata[:,12]
e_pmDEC=rdata[:,13]
Tpmllbb=bovy_coords.pmrapmdec_to_pmllpmbb(e_pmRA,e_pmDEC \
   ,RA_true,DEC_true,degree=True,epoch=2000.0)
e_pmGLON=Tpmllbb[:,0]
e_pmGLAT=Tpmllbb[:,1]

tbhdu = pyfits.BinTableHDU.from_columns([\
  pyfits.Column(name='RA_true',unit='(degree)',format='D',array=rdata[:,0]),\
  pyfits.Column(name='DEC_true',unit='(degree)',format='D',array=rdata[:,1]),\
  pyfits.Column(name='Plx_true',unit='(mas)',format='D',array=rdata[:,2]),\
  pyfits.Column(name='pmRA_true',unit='(mas/yr)',format='D',array=rdata[:,3]),\
  pyfits.Column(name='pmDEC_true',unit='(mas/yr)',format='D',array=rdata[:,4]),\
  pyfits.Column(name='HRV_true',unit='(km/s)',format='D',array=rdata[:,5]),\
# observed
  pyfits.Column(name='RA_obs',unit='(degree)',format='D',array=rdata[:,6]),\
  pyfits.Column(name='DEC_obs',unit='(degree)',format='D',array=rdata[:,7]),\
  pyfits.Column(name='Plx_obs',unit='(mas)',format='D',array=rdata[:,8]),\
  pyfits.Column(name='pmRA_obs',unit='(mas/yr)',format='D',array=rdata[:,9]),\
  pyfits.Column(name='pmDEC_obs',unit='(mas/yr)',format='D',array=rdata[:,10]),\
  pyfits.Column(name='HRV_obs',unit='(km/s)',format='D',array=rdata[:,11]),\
# error
  pyfits.Column(name='e_RA',unit='(degree)',format='D',array=rdata[:,12]),\
  pyfits.Column(name='e_DEC',unit='(degree)',format='D',array=rdata[:,13]),\
  pyfits.Column(name='e_Plx',unit='(mas)',format='D',array=rdata[:,14]),\
  pyfits.Column(name='e_pmRA',unit='(mas/yr)',format='D',array=rdata[:,15]),\
  pyfits.Column(name='e_pmDEC',unit='(mas/yr)',format='D',array=rdata[:,16]),\
  pyfits.Column(name='e_HRV',unit='(km/s)',format='D',array=rdata[:,17]),\
# True
  pyfits.Column(name='G_true',unit='(mag)',format='D',array=rdata[:,18]),\
  pyfits.Column(name='G_BP_true',unit='(mag)',format='D',array=rdata[:,20]),\
  pyfits.Column(name='G_RP_true',unit='(mag)',format='D',array=rdata[:,21]),\
# Observed
  pyfits.Column(name='G_obs',unit='(mag)',format='D',array=rdata[:,22]),\
  pyfits.Column(name='G_BP_obs',unit='(mag)',format='D',array=rdata[:,24]),\
  pyfits.Column(name='G_RP_obs',unit='(mag)',format='D',array=rdata[:,25]),\
# Error
  pyfits.Column(name='e_G',unit='(mag)',format='D',array=rdata[:,26]),\
  pyfits.Column(name='e_G_BP',unit='(mag)',format='D',array=rdata[:,28]),\
  pyfits.Column(name='e_G_RP',unit='(mag)',format='D',array=rdata[:,29]),\
# True
  pyfits.Column(name='Teff_true',unit='(K)',format='D',array=rdata[:,30]),\
  pyfits.Column(name='logg_true',unit='(dex)',format='D',array=rdata[:,31]),\
  pyfits.Column(name='[Fe/H]_true',unit='(dex)',format='D',array=rdata[:,32]),\
  pyfits.Column(name='Av_true',unit='(mag)',format='D',array=rdata[:,33]),\
# Observed
  pyfits.Column(name='Teff_obs',unit='(K)',format='D',array=rdata[:,34]),\
  pyfits.Column(name='logg_obs',unit='(dex)',format='D',array=rdata[:,35]),\
  pyfits.Column(name='[Fe/H]_obs',unit='(dex)',format='D',array=rdata[:,36]),\
  pyfits.Column(name='Av_obs',unit='(mag)',format='D',array=rdata[:,37]),\
# Errors
  pyfits.Column(name='e_Teff',unit='(K)',format='D',array=rdata[:,38]),\
  pyfits.Column(name='e_logg',unit='(dex)',format='D',array=rdata[:,39]),\
  pyfits.Column(name='e_[Fe/H]',unit='(dex)',format='D',array=rdata[:,40]),\
  pyfits.Column(name='e_Av',unit='(mag)',format='D',array=rdata[:,41]),\
# V,VI,GRVS
  pyfits.Column(name='V',unit='(mag)',format='D',array=rdata[:,42]),\
  pyfits.Column(name='V-I',unit='(mag)',format='D',array=rdata[:,43]),\
  pyfits.Column(name='G_RVS',unit='(mag)',format='D',array=rdata[:,44]),\
# GLON and GLAT only true
  pyfits.Column(name='GLON_true',unit='(degree)',format='D',array=GLON_true),\
  pyfits.Column(name='GLAT_true',unit='(degree)',format='D',array=GLAT_true),\
# pmGLON and pmGLAT true
  pyfits.Column(name='pmGLON_true',unit='(mas/yr)',format='D',\
                array=pmGLON_true),\
  pyfits.Column(name='pmGLAT_true',unit='(mas/yr)',format='D',\
                array=pmGLAT_true),\
# observed
  pyfits.Column(name='pmGLON_obs',unit='(mas/yr)',format='D',\
                array=pmGLON_obs),\
  pyfits.Column(name='pmGLAT_obs',unit='(mas/yr)',format='D',\
                array=pmGLAT_obs),\
# Errors
  pyfits.Column(name='e_pmGLON',unit='(mas/yr)',format='D',\
    array=e_pmGLON),\
  pyfits.Column(name='e_pmGLAT',unit='(mas/yr)',format='D',\
    array=e_pmGLAT)])
tbhdu.writeto(outfile,clobber=True)

