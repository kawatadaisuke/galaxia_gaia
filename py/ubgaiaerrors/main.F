c---------------------------------------------------------------------
c main.f 
c   modified from galactocentric_observed.f
c       M. Romero-Gomez, F. Figueras, T. Antoja
c       ICCUB-IEEC
c       October 2012
c
c       input: from pyana/plotrc.py
c
c       output: xf vector with galactocentric coordinates affected
c               with Gaia errors, and written in the output file
c
c       The code assumes the Sun is located on the x-negative axis, its
c       (xsun,0.0,zsun) xsun>0
c       vlsr: LSR rotation velocity
c       and moving with respect to LSR (usun,vsun,wsun)
c
c       The photometry of the population is specified in const_pop.h
c
c---------------------------------------------------------------------
      implicit none

      integer np
      integer i,j,ip,jflag,idr
      double precision month,CAfactor
      double precision x,y,z,vx,vy,vz,mp,mag,G,GRVS,GBP,GRP
      double precision alpha,delta,uh,vh,wh,dist,xpi,gl,gb
c position and velocity of Sun
      double precision xsun,zsun,usun,vsun,wsun,vlsr
c luminosity and colour
      double precision v,vi,Av
c stellar parameter
      double precision Teff,FeH,xlogg,lrp,ri,ro,Age
c random seed
      integer idum,ndr    
c observable
      double precision a(6),ao(6),ae(6),p(4),po(4),pe(4),ap(4)
     & ,apo(4),ape(4)
      double precision mua,mud,vrad
      double precision pi
      PARAMETER (pi = 3.14159265358979324d0)

c *** for RC stars ***
c      mag=1.2714d0
c      VI=1.0d0
c      Teff=3850.0d0

c *** mission length in months, 60=5 years ***
c *** 22 months for Gaia DR2 (testing ongoing) ***
      month=22.0d0
c *** Calibration astrometric factor:
c *** CAfactor=1.0d0 means no calibration errors added ***
c *** CAfactor=1.2d0 gives 20% calibration error (on the safe side) ***
c      CAfactor=1.2d0
c DR2 seems to be better than expected.
      CAfactor=1.0d0
c Introduce Gaia errors
c        jflag=-1                          !weighted errors
       jflag= 1                          !mean errors after commissioning, but DR2 seems to be better than this

c set random seed
      idum=-1111
      
c *** input file
      open(50,file='../ubgaiaein.bin',status='old'
     &  ,form='unformatted')
      read(50) np
      write(6,*) 'Number of Stars=',np

c *** output file
      open(60,file='../gaiaei-out.bin',status='unknown'
     &  ,form='unformatted')
      write(60) np
      close(60)

      open(60,file='../gaiaed-out.bin',status='unknown'
     &  ,form='unformatted')

c      open(70,file='./test.dat',status='unknown')
      do j=0,np-1
        read(50) alpha,delta,dist,mua,mud,vrad,Teff,xlogg,FeH,Av
     &     ,vi,v,Age
c        write(70,'(2(1pE13.5))') Teff,xlogg
c distance (pc) -> parallax (mas)
        xpi=1000.0d0/dist
c         a(1): Equatorial heliocentric right ascension (units: radians)
        a(1)=alpha*pi/180.0d0
c         a(2): Equatorial heliocentric declination: delta (units: radians) 
        a(2)=delta*pi/180.0d0 
c         a(3): Parallax: pi (units: mas)
        a(3)=xpi
c         a(4): Equatorial proper motions in right ascension in true arcs on the sky: mua_{*}=mua*cos(delta) (units: mas/yr)
        a(4)=mua*dcos(delta)
c         a(5): Equatorial proper motions: mud (units: mas/yr)
        a(5)=mud
c         a(6): Radial velocity: Vr (units: km/s) 
        a(6)=vrad

c Assign atmospheric parameters to the source
        ap(1)=Teff    
        ap(2)=xlogg  
        ap(3)=FeH
        ap(4)=Av


        call Gaia_errors(month,CAfactor,jflag,V,VI,a,ao,ae,p,po
     &    ,pe,ap,apo,ape,idum)

c truearc value -> normal value
        if(dcos(delta).eq.0.0d0) then
          ae(1)=ae(1)/dcos(delta)    
          a(4)=a(4)/dcos(delta)    
          ae(4)=ae(4)/dcos(delta)    
          ao(4)=ao(4)/dcos(delta)    
        endif

c Table 3 of Jordi et al. (2010, A&A, 523, 48
        G=V-0.0257-0.0924d0*VI-0.1623d0*(VI**2)+0.0090d0*(VI**3)
        GRVS=V-0.0119d0-1.2092d0*VI+0.0188d0*(VI**2)+0.0005d0*(VI**3)
        GBP=V-(0.0643d0-0.3266d0*VI+0.0887d0*(VI**2)-0.0050d0*(VI**3))
        GRP=V-(-0.0017d0+0.8794d0*VI+0.0273d0*(VI**2)-0.0008d0*(VI**3))
        write(60) (a(i),i=1,6),(ao(i),i=1,6),(ae(i),i=1,6)
     &     ,(p(i),i=1,4),(po(i),i=1,4),(pe(i),i=1,4)
     &     ,(ap(i),i=1,4),(apo(i),i=1,4),(ape(i),i=1,4)
     &     ,V,VI,GRVS,Age
      enddo
      close(60)
c      close(70)

      stop
      end

