c ---------------------------------------------------------------------
c
c Subroutine that computes and assign Gaia errors to observables
c Errors are computed following Gaia Science Performance web page:
c http://www.rssd.esa.int/index.php?project=GAIA&page=Science_Performance
c
c Universitat de Barcelona
c Contact: Merce Romero-Gomez, Josep Manel Carrasco 
c
c ---------------------------------------------------------------------

c Modified by D. Kawata (MSSL, UCL): 24 January 2018
c added BP RP photmetry output

c Inputs:
c ---------

c   * External tables:
c     ----------------
c	'gfactor-Oct2012.dat': geometrical factors
c 	'TableVr-Oct2012.dat': Radial Velocity coefficients (pre launch)
c       'TableVr-Oct2014.dat': Update to post launch errors

c	Input files for the definition of constants: 
c   		const_math.h : mathematical constants
c   		const_ast.h  : astronomical constants 


c   * Parameters from the star:
c     -------------------------

c     jflag:

c     	The code offers two options depending of the type of errors 
c     	you want to apply (see website): 
c     		jflag=1, errors computed from mean geometrical factors (see website)
c     		jflag=-1,errors computed considering the scanning law of the satellite 
c             		(the geometrical factors are then computed from the ecliptic
c             		latitude and the number transits)

c     V: V apparent magnitude 
c     VI: observed (V-I) colour index (affected by extinction), in magnitudes

c     a(i): Astrometric actual values for the star (not affected by errors)
c         a(1): Equatorial heliocentric right ascension (units: radians)
c         a(2): Equatorial heliocentric declination: delta (units: radians) 
c         a(3): Parallax: pi (units: mas)
c         a(4): Equatorial proper motions in right ascension in true arcs on the sky: mua_{*}=mua*cos(delta) (units: mas/yr)
c         a(5): Equatorial proper motions: mud (units: mas/yr)
c         a(6): Radial velocity: Vr (units: km/s) 

c     ap(i): Atmospheric parameters actual values for the star (no affected by errors) 
c         ap(1): Effective temperature (units: kelvin)
c         ap(2): surface gravitye (logg) (units: dex) 
c         ap(3): [Fe/H] (units: dex) 
c         ap(4): Absorption A0 (Av~A0, see Bailer-Jones et al., 2011)



c Output values: 
c -------------
c Parameters of the star as observed by Gaia (affected by errors)
c Gaia errors assigned to each parameter:  

c The output parameters of the star as observed by Gaia are: 

c   Astrometric (+Vr) data affected by errors  
c     ao(i): Astrometric values for the star affected by Gaia errors
c         ao(1): Observed Equatorial heliocentric right ascension (units: radians)
c         ao(2): Observed Equatorial heliocentric declination: delta (units: radians) 
c         ao(3): Observed Parallax: pi (units: mas)
c         ao(4): Observed Equatorial proper motions in right ascension in true arcs on the sky: mua_{*}=mua*cos(delta) (units: mas/yr)
c         ao(5): Observed Equatorial proper motions: mud (units: mas/yr)
c         ao(6): Observed Radial velocity: Vr (units: km/s) 
c     ae(i): Gaia standard deviation (1-sigma) 
c         ae(1): Standard deviation in Equatorial heliocentric right ascension in true arcs on the sky:: alpha_{*}=alpha*cos(delta)  (units: mas)
c         ae(2): Standard deviation in Equatorial heliocentric declination: delta (units: mas) 
c         ae(3): Standard deviation in Parallax: pi (units: mas)
c         ae(4): Standard deviation in Equatorial proper motions in right ascension in true arcs on the sky: mua_{*}=mua*cos(delta) (units: mas/yr)
c         ae(5): Standard deviation in Equatorial proper motions: mud (units: mas/yr)
c         ae(6): Standard deviation in Radial velocity: Vr (units: km/s) 

c   Photometric data affected by errors (G, G_BP, G_RP, G_RVS)
c     p(i): Photometric actual values for the star (no affected by errors)
c         p(1): G magnitude (units: magnitudes) 
c         p(2): G_BP-G_RP   (units: magnitudes) 
c DK added G_BP and G_RP output
c         p(3): G_BP (units: magnitudes) 
c         p(4): G_RP   (units: magnitudes) 
c     po(i): Photometric values of the star affected by Gaia errors
c         po(1): Observed G magnitude (units: magnitudes) 
c         po(2): Observed G_BP-G_RP   (units: magnitudes) 
c         po(3): Observed G_BP (units: magnitudes) 
c         po(4): Observed G_RP   (units: magnitudes) 
c     pe(i): Gaia standard deviation (errors, 1-sigma)
c         pe(1): Standard deviation in G magnitude (units: magnitudes) 
c         pe(2): Standard deviation in G_BP-G_RP   (units: magnitudes) 
c         pe(3): Standard deviation in G_BP (units: magnitudes) 
c         pe(4): Standard deviation in G_RP   (units: magnitudes) 

c   Atmospheric parameters data affected by errors: 
c     apo(i): Atmospheric parameters of the star affected by Gaia errors
c         apo(1): Observed Effective temperature (units: kelvin)
c         apo(2): Observed surface gravity (logg) (units: dex) 
c         apo(3): Observed [Fe/H] (units: dex) 
c         apo(4): Observed Absorption A0 (Av~A0, see Bailer-Jones et al., 2011) (units: magnitudes) 
c     ape(i): Gaia standard deviation (1-sigma)
c         ap(1): Standard deviation in Effective temperature (units: kelvin)
c         ap(2): Standard deviation in surface gravity (logg) (units: dex) 
c         ap(3): Standard deviation in [Fe/H] (units: dex) 
c         ap(4): Standard deviation in Absorption A0 (Av~A0, see Bailer-Jones et al., 2011) (units: magnitudes) 



c **********************************************************************
        subroutine Gaia_errors(month,CAfactor,jflag,V,VI,a,ao,ae,p,po,pe
     &             ,ap,apo,ape,Llavor)
        implicit real*8 (a-h,o-z)
        dimension a(6),ao(6),ae(6),p(4),po(4),pe(4),g(5)
        dimension ap(4),apo(4),ape(4)
        dimension xvi(11),xavr(11),xbvr(11)
c        dimension xsb(20),xNOBS(20),xga(20),xgd(20),xgpi(20),xgmua(20),
c     &            xgmud(20)
        dimension xsb(20),xga(20),xgd(20),xgpi(20),xgmua(20),
     &            xgmud(20)
c        implicit none
c        double precision a(6),ao(6),ae(6),p(4),po(4),pe(4),g(5),
c     &       ap(4),apo(4),ape(4),xvi(11),xavr(11),xbvr(11),
c     &       xsb(20),xNOBS(20),xga(20),xgd(20),xgpi(20),xgmua(20),
c     &       xgmud(20)
c        double precision V,VI,empi,gmeana,gmeand,gmeanmua,gmeanmud,
c     &        gmeanpi,sbeta
c        integer jflag,i,j

        double precision month,CAfactor
        real*8 gasdev
        character*10 xchar
        logical ifirst
        data ifirst/.true./
c        save xsb,xNOBS,xga,xgd,xgpi,xgmua,xgmud,xnmean,gmeana,gmeand,
c     &       gmeanpi,gmeanmua,gmeanmud
        save xsb,xga,xgd,xgpi,xgmua,xgmud,gmeana,gmeand,
     &       gmeanpi,gmeanmua,gmeanmud

        INCLUDE 'const_math.h'
        INCLUDE 'const_ast.h'
        COMMON / tableVr / xvi, xavr, xbvr 
        integer Llavor

        factorL=sqrt(60.d0/month)
        L=month/12.d0

c        write(6,*) factorL,L,CAfactor

        if(ifirst)then
          write(6,*) ' reading gfactor-Jun2013'
c        open(3,file='gfactor-Oct2012.dat')
        open(3,file='gfactor-Jun2013.dat')
            read(3,*)
            do i=1,20
              read(3,*)xsb(i),xga(i),xgd(i),xgpi(i),
     &                 xgmua(i),xgmud(i)

c              write(6,'(6f6.3)') xsb(i),xga(i),xgd(i),xgpi(i),
c     &                 xgmua(i),xgmud(i)

            enddo
            read(3,*) xchar,gmeana,gmeand,gmeanpi,gmeanmua,
     &               gmeanmud

c              write(6,'(a4,5f6.3)') 'mean',gmeana,gmeand,gmeanpi
c     &              ,gmeanmua,gmeanmud

            close(3)

       open(3,file='TableVr-Jun2015.dat')
            read(3,*)
            do i=1,11
              read(3,*)xvi(i),xavr(i),xbvr(i)
            enddo
            close(3)
        ifirst=.false.
        endif

c Get Gaia errors in parallax,magnitudes and Vr 

        call ErrorsPiMag(factorL,CAfactor,v,vi,p,pe,empi)
        call ErrorsVr(v,vi,ae(6),GRVS)

c       Output errors on parallax and proper motion will be in mas
        empi = empi /1.0e3

c       Mean errors (from mean g-geometrical factor)
        if(jflag.gt.0)then
          ae(3)=empi
          ae(1)=gmeana*ae(3)
          ae(2)=gmeand*ae(3)
          ae(4)=gmeanmua*ae(3)*(factorL**2)
          ae(5)=gmeanmud*ae(3)*(factorL**2)
        endif

c       Errors depending on the scanning law 
c       (number of passages and the geometrical factor)
        if(jflag.lt.0)then
           
          sbeta=dabs(0.9175d0*dsin(a(2))
     &      -0.3978d0*dcos(a(2))*dsin(a(1)))

c          call lininter(20,xsb,xNobs,sbeta,gns)
          call lininter(20,xsb,xga,sbeta,g(1))
          call lininter(20,xsb,xgd,sbeta,g(2))
          call lininter(20,xsb,xgpi,sbeta,g(3))
          call lininter(20,xsb,xgmua,sbeta,g(4))
          call lininter(20,xsb,xgmud,sbeta,g(5))

c          xnob=dsqrt(xnmean/gns)

          do i=1,5
c          ae(i)=xnob*g(i)*empi
            ae(i)=g(i)*empi
c            write(6,*) ae(4),ae(5)
          enddo 
c DK added parallax and proper motion error dependence of month
c Anthony Brown Lund talk http://www.astro.lu.se/gaia2017/slides/Brown.pdf
c parallax improve t^-0.5
c proper motion improve t^-1.5
          ae(3)=ae(3)*factorL
          ae(4)=ae(4)*(factorL**3)
          ae(5)=ae(5)*(factorL**3)

        endif


c       Computation of the errors in the atmospheric parameters
c       See Liu et al., 2012
c       A second order polinomial on G has been fitted to the
c       |Teff_eaneas_pq - Teff_real|

        call errorsAp(p(1),ape)
	
c       Computation of the observed astrometric quantities 

c       The error in right ascension ae(1) denotes true arc on the sky
c       so the right ascension shall be converted to that before the 
c       random error is assigned 
c       alpha_{*}=alpha*cos(delta) 
 
        a(1)=a(1)*dcos(a(2))

c	The errors in alpha* and  delta are in mas whereas alpha and delta are
c       in radians, so (alpha*,delta) are converted from radians to mas  
c       before the random error is assigned 
 
	a(1)=a(1)/mas
	a(2)=a(2)/mas
	
        do i=1,6
           ao(i)=gasdev(a(i),ae(i),Llavor)
c           write(6,*) i,a(i),ao(i),ae(i),Llavor
        enddo 
c        write(6,*) 'change'
c       The true and observed (alpha*,delta) are converted from mas to radians 
	a(1)=a(1)*mas
	a(2)=a(2)*mas
	ao(1)=ao(1)*mas
	ao(2)=ao(2)*mas

c        a(4)=a(4)*mas
c        ao(4)=ao(4)*mas
c        ae(4)=ae(4)*mas
c       and the alpha_{*} is converted to alpha
        ao(1)=ao(1)/dcos(ao(2))
        a(1)=a(1)/dcos(a(2))
        ae(1)=a(1)/dcos(a(2))

c        a(4)=a(4)/dcos(a(2))

c        ao(4)=ao(4)/dcos(a(2))
c        ae(4)=ae(4)/dcos(a(2))
c        a(4)=a(4)/mas
c        ao(4)=ao(4)/mas
c        ae(4)=ae(4)/mas
c       Computation of the observed photometric quantities 
        do i=1,4
        po(i)=gasdev(p(i),pe(i),Llavor)
        enddo

c       Computation of the observed atmospheric parameters
        do i=1,4
        apo(i)=gasdev(ap(i),ape(i),Llavor)
        enddo

        return
	end


c******************************************************************
c     Subroutine ErrorsAp
c
c     A second order polinomial on G has been fitted to the
c     |Teff_eaneas_pq - Teff_real|
c     Commented, there is the option to fitting a second order
c     polinomial to the sigmas of the Aeneas pq-model 
c     (see tableA4 of Liu & Bailer-Jones 2012)
c    
c     input data
c     p(1)=G           Gaia apparent magnitude (G)

c     output data
c     ape(1):          Uncertainty in Teff (Kelvin) 
c     ape(2):          Uncertainty in logg (dex) 
c     ape(3):          Uncertainty in [Fe/H] (dex) 
c     ape(4):          Uncertainty in A0 (mag) 

c------------------------------------------------------------------

      subroutine errorsAp(G,ape)

      implicit real*8 (a-h,o-z)
      dimension ape(4)


c     Fitting |Teff_eaneas_pq - Teff_real|
      ape(1)= 630.3d0   - 107.54d0  *G + 4.9682d0    *G**2
      ape(2)= 0.49865d0 - 0.04432d0 *G + 0.0017055d0 *G**2
      ape(3)= 0.81192d0 - 0.12226d0 *G + 0.0056669d0 *G**2
      ape(4)= 0.56387d0 - 0.093435d0*G + 0.0042024d0 *G**2

c     Fitting sigTeff (Eaneas pq-model)
c      ape(1)=67.289d0   - 7.8492d0  *G + 0.34687d0   *G**2
c      ape(2)= 0.31915d0 - 0.044319d0*G + 0.0020501d0 *G**2
c      ape(3)= 0.71709d0 - 0.115d0   *G + 0.0050649d0 *G**2
c      ape(4)= 0.17787d0 - 0.026371d0*G + 0.0011578d0 *G**2

      return 
      end


c******************************************************************
c     Subroutine ErrorsPiMag
c
c     input data
c     v=V           Johnson V apparent magnitude
c     vi=(V-I)c     Johnson-Cousins V-I colour

c     output data
c     p(1)=G                true Gaia G apparent magnitude
c     p(2)=G_BP-G_RP        true Gaia colour
c     pe(1)=sigG            uncertainty in G (in mag units)
c     pe(2)=sigBPRP         uncertainty in G_BP-G_RP colour (in mag units)
c     ae(3)=empi            uncertainty in parallax (in microarcsec units)

c------------------------------------------------------------------
      subroutine ErrorsPiMag(factorL,CAfactor,v,vi,p,pe,empi)

      implicit real*8 (a-h,o-z)
      dimension p(4),pe(4)


c---- Gaia photometry,

      p(2)=-0.0660d0+1.2061d0*vi-0.0614d0*vi*vi+0.0041d0*vi**3
      p(1)=V-0.0257d0-0.0924d0*vi-0.1623d0*vi*vi+0.0090d0*vi**3 
      if(p(1).gt.20.) then
         empi=999999.999
         pe(1)=99.999
         pe(2)=99.999
         return
      endif

c      factorL=1.0d0
c      CAfactor=1.0d0
c Astrometric Error
c After Comissioning Nov 2014
      z=MAX(10**(0.4*(12.0-15.0)), 10**(0.4*(p(1)-15.0)))
      empi=dsqrt(-1.631+680.766*z+32.732*z**2)*(0.986+(1.0-0.986)*vi)

      empi=empi*factorL*CAfactor
      if(p(1).lt.15.0d0.and.empi.lt.40.0d0) then
        empi=40.0d0
      endif

c Photometric Errors
      pe(1)=0.001*dsqrt(0.04895*z**2+1.8633*z+0.0001985)
      pe(1)=(pe(1)/dsqrt(70.d0))*factorL

      z=MAX(10**(0.4*(11.0-15.0)), 10**(0.4*(p(1)-15.0)))

      a=-0.000562*vi**3+0.044390*vi**2+0.355123*vi+1.043270
      b=-0.000400*vi**3+0.018878*vi**2+0.195768*vi+1.465592
      c= 0.000262*vi**3+0.060769*vi**2-0.205807*vi-1.866968
      sigBP=0.001*dsqrt((10.**a)*(z**2)+(10.**b)*z+10.**c)
      sigBP=sigBP/dsqrt(70.d0) * factorL
      pe(3)=sigBP
      
      a=-0.007597*vi**3+0.114126*vi**2-0.636628*vi+1.615927
      b=-0.003803*vi**3+0.057112*vi**2-0.318499*vi+1.783906
      c=-0.001923*vi**3+0.027352*vi**2-0.091569*vi-3.042268
      sigRP=0.001*dsqrt((10.**a)*(z**2)+(10.**b)*z+10.**c)
      sigRP=sigRP/dsqrt(70.d0)* factorL
      pe(4)=sigRP

c     We assume there is no correlation among colours
      pe(2)=dsqrt(sigBP**2+sigRP**2) 

      return
      end



c******************************************************************
c     Subroutine ErrorsVr
c
c     input data
c     v=V           Johnson V apparent magnitude
c     vi=(V-I)c     Johnson-Cousins V-I colour

c     output data
c     ae(6)=Sigr            uncertainty in Vr (in km/s units)
c------------------------------------------------------------------
      subroutine ErrorsVr(v,vi,sigVr,GRVS)

      implicit real*8 (a-h,o-z)
      dimension xavr(11),xbvr(11),xvi(11)
      COMMON / tableVr / xvi, xavr, xbvr 
c
c     Vr errors 
c     
     
      GRVS=V-0.0119D0-1.2092D0*VI+0.0188D0*VI*VI+0.0005D0*VI*VI*VI
      call lininter(11,xvi,xavr,vi,avr)
      call lininter(11,xvi,xbvr,vi,bvr)

      sigVr=1.d0+bvr*exp(avr*(V-12.7d0))
      return
      end

c **********************************************************************
c          subroutine lininter(n,xa,ya,x,y)
c input:
c   xa(n) and ya(n) such that ya(i)=f(xa(i)), and xa ordered in 
c                   increasing sense
c   x, point where you want to know the corresponding y=f(x)
c output:
c   y
C  --------------------------------------------------------------------

           subroutine lininter(n,xa,ya,x,y)
           implicit real*8(a-h,o-z)
           dimension xa(n),ya(n)
          
           if(x.gt.xa(n))then
             y=ya(n)
             return
           else if(x.lt.xa(1))then
             y=ya(1)
             return
           else if(xa(1).le.x.and.x.le.xa(n))then
             i=1
             do while (xa(i).lt.x)
               i=i+1
             enddo
             i=i-1
             if(i.le.0) then
               write(6,*) ' Error in linter: i=0'
               write(6,*) ' n,x=',n,x
             else if(i.ge.n) then
               write(6,*) ' Error in linter: i=n'
               write(6,*) ' n,x=',n,x
             endif
             x0=xa(i)
             x1=xa(i+1)

c             write(6,*) ' x,x0,x1=',x,x0,x1

             y0=ya(i)
             y1=ya(i+1)
             y=((x-x0)/(x1-x0))*y1+((x1-x)/(x1-x0))*y0
             return
           endif

           end

c******************************************************************
C  Function ran1 
C  --------------------------------------------------------------------

      FUNCTION ran1(idum)
      INTEGER idum, IA, IM, IQ, IR, NTAB, NDIV
      REAL ran1, AM, EPS, RNMX
      PARAMETER (IA=16807, IM= 2147483647, AM=1./IM,IQ=127773,IR=2836,
     *     NTAB=32,NDIV=1+(IM-1)/NTAB, EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
          idum=max(-idum,1)
          do j=NTAB+8,1,-1
             k=idum/IQ
             idum=IA*(idum-k*IQ)-IR*k
             if (idum.lt.0) idum=idum+IM
             if (j.le.NTAB) iv(j)=idum
          enddo
          iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      end



c******************************************************************

c     subroutine get_fi
c input: (alfa,delta) in radians
c output: fi parallactic angle in radians

C  --------------------------------------------------------------------
      subroutine get_fi(alfa,delta,fi)
      implicit real*8(a-h,o-z)
      include 'const_math.h'
      include 'const_ast.h'

      x=Dcos(Dg)*Dsin(Delta)*Dsin(Alfa-Ag)+Dsin(Dg)*Dcos(Delta)
      y=Dcos(Dg)*Dcos(Alfa-Ag)
c Parallactic angle
      Fi=Datan2(y,x)
      return
      end

c******************************************************************
c  Generation of a random deviate drawn from a Gaussian Distribution 
C  --------------------------------------------------------------------

      Function gasdev(X0,SigX,idum)
      Integer idum
      double precision gasdev,X0,SigX
      INTEGER iset
      double precision fac,gset,rsq,v1,v2
      real ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
 1        v1=2.d0*dble(ran1(idum))-1.d0
          v2=2.d0*dble(ran1(idum))-1.d0
          rsq=v1**2+v2**2
          if(rsq.ge.1.d0.or.rsq.eq.0.d0) goto 1
          fac=dsqrt(-2.d0*log(rsq)/rsq)
          gset=v1*fac
          gasdev=v2*fac*SigX + X0
c          write(6,*) 'GD1',v2,fac,SigX,X0,gasdev
          iset=1
       else
          gasdev=gset*SigX + X0
c          write(6,*) 'GD2',v2,fac,SigX,X0,gasdev
          iset=0
       endif
       return
       end

C  --------------------------------------------------------------------
c     SUBROUTINE carte_to_equatorial(X,Y,Z,xsun,zsun,ALPHA,DELTA)
c     Input cartesian coordinates : (X,Y,Z) in kpc
c           location of sun : (xsun,zsun) in kpc, xsun<0
c     Output: equatorial coordinates: (alpha,delta) in radiants
c     
c---------------------------------------------------------------------
      SUBROUTINE carte_to_equatorial(x,y,z,xsun,zsun,alpha,delta,gl,gb)
      implicit none
      include 'const_math.h'
      include 'const_ast.h'
      double precision x,y,z,xsun,zsun
      double precision alpha,delta
      double precision sx,sy,sz,gr,gb,gl,cb,sb,cl,sl,sd,yk,xk,gr2
c      parameter (pi=3.14159265358979324d0)

c XYZ, centered on the Sun
      sx=x-xsun
      sy=y
      sz=z-zsun

c R,l,b heliocentric
       gr=dsqrt(sx*sx+sy*sy+sz*sz)
       gb=asin(sz/gr)
       gr2=dsqrt(sx**2+sy**2)
      if(sy.gt.0.0d0) then
        gl=dacos(sx/gr2) 
      else
        gl=2.0d0*pi-dacos(sx/gr2) 
      endif

c alpha, delta heliocentric
      cb=Dcos(gb)
      sb=Dsin(gb)
      cl=Dcos(gl)
      sl=Dsin(gl)
      
      sd=cb*Dsin(gl-a2)*Dsin(a1)+
     +     sb*Dcos(a1)
      Delta=Dasin(sd)
      
      yk=cb*Dsin(gl-a2)*Dcos(a1)-
     .     sb*Dsin(a1)
      xk=cb*Dcos(gl-a2)
      Alpha=mod(Datan2(yk,xk) +Ag,2.d0*pi)

      RETURN
      END


      SUBROUTINE UVWH_to_equatorials(U,V,W,A)
C*************************************************************************
C USING TRIGONOMETRIC MATRIX COMPUTE (VR,MUA,MUD) FROM (U,V,W)HELIOCENTRIC
c input velocities in km/s, parallax, xpi, in arcsec, angles in radians
C a(1)=alpha
C a(2)=delta    the first three are input variables
C a(3)=xpi    
C a(4)=mua
C a(5)=mud      the last three are output variables
C a(6)=Vr
c output proper motions in "/yr and radial velocity in km/s
C*************************************************************************
      implicit none
      INCLUDE 'const_math.h'
      INCLUDE 'const_ast.h'
      double precision a(6)
      double precision u,v,w

      double precision fi,sfi,cfi,gl,gb,cb,sb,cl,sl,gr
      double precision b11,b12,b13,b21,b22,b23,b31,b32,b33

      call get_fi(a(1),a(2),fi)
      SFI=DSIN(FI)
      CFI=DCOS(FI)
      call equ_to_gal(a(1),a(2),gl,gb)
C APPLY TRANSFORMATION MATRIX TO (U,V,W)INPUT
      cb=Dcos(gb)
      sb=Dsin(gb)
      cl=Dcos(gl)
      sl=Dsin(gl)

      B11=CL*CB 
      B12=CL*SB*SFI-SL*CFI
      B13=-(CL*SB*CFI+SL*SFI)
      B21=SL*CB
      B22=SL*SB*SFI+CL*CFI
      B23=-(SL*SB*CFI-CL*SFI)
      B31=SB
      B32=-CB*SFI
      B33=CB*CFI

      
      a(6)=B11*U+B21*V+B31*W
      a(4)=B12*U+B22*V+B32*W
      a(5)=B13*U+B23*V+B33*W
      gr=1.d0/a(3)
C CORRECT UNITS 
      a(4)=a(4)/(kt*gr)
      a(5)=a(5)/(kt*gr)

      RETURN
      END


c **********************************************************************

c subroutine equatorial_to_galactic(alfa,delta,gl,gb)
c input: equatorial coordinates (alfa,delta)
c output: galactic coordinates (xl,xb)
c Units: all angles l,b,alpha,delta in radians
c -------------------------------------------------------------

      subroutine equ_to_gal(alfa,delta,gl,gb)
      implicit none
      include 'const_math.h'
      include 'const_ast.h'
      double precision alfa,delta,gl,gb
      double precision t6,t36,rd,lq,alfag,decg,sdg,cdg,sd,cd,sb,cb
     & ,slpl,clpl,lpl
c      parameter (pi=3.14159265358979324d0)

      T6=60.d0
      T36=3600.d0
      RD=PI/180.d0
      LQ=122.92851d0*RD
      ALFAG=(12.d0+51.d0/T6+26.2755d0/T36)*15.d0*RD
      DECG=(27.d0+7.d0/T6+41.704d0/T36)*RD
      SDG=dSIN(DECG)
      CDG=dCOS(DECG)
 
      sd=dsin(delta)
      cd=dcos(delta)

      sb=sdg*sd+cdg*cd*dcos(alfa-alfag)
      cb=sqrt(1.d0-sb**2.d0)
      gb=datan2(sb,cb)
      slpl=(cd*dsin(alfa-alfag))/cb
      clpl=(cdg*sd-sdg*cd*dcos(alfa-alfag))/cb
      lpl=datan2(slpl,clpl)
      if(lpl.le.0.d0)lpl=lpl+2.d0*pi
      gl=lq-lpl

      return
      end

c ****************************************************************
c Subrountine coord 
c ****************************************************************
      subroutine coord(ao,ae,xo_p,yo_p,zo_p,Uo,Vo,Wo,
     &     Ue,Ve,We,lo_p,bo_p,vlo_p,vbo_p)

      implicit none
      double precision SI_tmateqgal(0:2,0:2),SI_tmatgaleq(0:2,0:2)
      double precision d,xeq,yeq,zeq,rxy,xo_p,yo_p,zo_p
     &     ,vxo_p,vyo_p,vzo_p,Uo,Vo,Wo,A(3),r(3),ao(6),b(6),vl_p
     &     ,vb_p,vlo_p,vbo_p,l1,l2,l3,c1,c2,c3,vl2,vlo2
     &     ,AA(0:2,0:2),BB(0:2,0:2),CC(0:2,0:2),r1,r2,r3
     &     ,s1,s2,s3,sfac,Ue,Ve,We,ae(6),vrots
     &     ,phi,sinphi,lo_p,bo_p,Tao4,Tao5
      INCLUDE 'const_math.h'

c galactic (Cartesian, Solar origin) -> heliocentlic (Cartesian) 
      SI_tmateqgal(0,0)=-0.054875290d0
      SI_tmateqgal(1,0)=0.494109674d0
      SI_tmateqgal(2,0)=-0.867666026d0
      SI_tmateqgal(0,1)=-0.873437003d0
      SI_tmateqgal(1,1)=-0.444830323d0
      SI_tmateqgal(2,1)=-0.198075198d0
      SI_tmateqgal(0,2)=-0.483834719d0
      SI_tmateqgal(1,2)=0.746981441d0
      SI_tmateqgal(2,2)=0.455985407d0
c heliocentlic (Cartesian) -> galactic (Cartesian, Solar origin) 
      SI_tmatgaleq(0,0)=-0.054877639d0
      SI_tmatgaleq(1,0)=-0.873437219d0
      SI_tmatgaleq(2,0)=-0.483835031d0
      SI_tmatgaleq(0,1)=0.494110185d0
      SI_tmatgaleq(1,1)=-0.444829280d0
      SI_tmatgaleq(2,1)=0.746982180d0
      SI_tmatgaleq(0,2)=-0.867665586d0
      SI_tmatgaleq(1,2)=-0.198076592d0
      SI_tmatgaleq(2,2)=0.455983866d0
c *** Calculate catesian error added positions ***
      d=1.0d0/(ao(3)*1.0e5)
      rxy=d*dcos(ao(2))
      xeq=rxy*dcos(ao(1))
      yeq=rxy*dsin(ao(1))
      zeq=d*dsin(ao(2))
      xo_p=SI_tmatgaleq(0,0)*xeq+SI_tmatgaleq(1,0)*yeq
     &     +SI_tmatgaleq(2,0)*zeq
      yo_p=SI_tmatgaleq(0,1)*xeq+SI_tmatgaleq(1,1)*yeq
     &     +SI_tmatgaleq(2,1)*zeq
      zo_p=SI_tmatgaleq(0,2)*xeq+SI_tmatgaleq(1,2)*yeq
     &     +SI_tmatgaleq(2,2)*zeq
c *** Cartesian velocities ***
      r(1)=ao(6)
      r(2)=ao(4)*4.74d0/ao(3)
      r(3)=ao(5)*4.74d0/ao(3)
      A(1)=r(1)*dcos(ao(1))*dcos(ao(2))-r(2)*dsin(ao(1))
     &     -r(3)*dcos(ao(1))*dsin(ao(2))
      A(2)=r(1)*dsin(ao(1))*dcos(ao(2))+r(2)*dcos(ao(1))
     &     -r(3)*dsin(ao(1))*dsin(ao(2))
      A(3)=r(1)*dsin(ao(2))+r(3)*dcos(ao(2))
      Uo=SI_tmatgaleq(0,0)*A(1)+SI_tmatgaleq(1,0)*A(2)+
     &     SI_tmatgaleq(2,0)*A(3)
      Vo=SI_tmatgaleq(0,1)*A(1)+SI_tmatgaleq(1,1)*A(2)+
     &     SI_tmatgaleq(2,1)*A(3)
      Wo=SI_tmatgaleq(0,2)*A(1)+SI_tmatgaleq(1,2)*A(2)+
     &     SI_tmatgaleq(2,2)*A(3)
c *** Cartesian velocity errors? (BETA - testing required) ***
      r1=ae(6)**2
      r2=((4.74d0/ao(3))**2)*(ae(4)**2+(ao(4)*ae(3)/ao(3))**2)
      r3=((4.74d0/ao(3))**2)*(ae(5)**2+(ao(5)*ae(3)/ao(3))**2)
      sfac=2.0d0*ao(4)*ao(5)*(4.74d0**2)*(ae(3)**2)/(ao(3)**4)
      AA(0,0)=dcos(ao(1))*dcos(ao(2))
      AA(0,1)=-1.0d0*dsin(ao(1))
      AA(0,2)=-1.0d0*dcos(ao(1))*dsin(ao(2))
      AA(1,0)=dsin(ao(1))*dcos(ao(2))
      AA(1,1)=dcos(ao(1))
      AA(1,2)=-1.0d0*dsin(ao(1))*dsin(ao(2))
      AA(2,0)=dsin(ao(2))
      AA(2,1)=0.0d0
      AA(2,2)=dcos(ao(2))
      BB(0,0)=SI_tmateqgal(0,0)*AA(0,0)+SI_tmateqgal(0,1)*AA(1,0)
     &     +SI_tmateqgal(0,2)*AA(2,0)
      BB(0,1)=SI_tmateqgal(0,0)*AA(0,1)+SI_tmateqgal(0,1)*AA(1,1)
     &     +SI_tmateqgal(0,2)*AA(2,1)
      BB(0,2)=SI_tmateqgal(0,0)*AA(0,2)+SI_tmateqgal(0,1)*AA(1,2)
     &     +SI_tmateqgal(0,2)*AA(2,2)
      BB(1,0)=SI_tmateqgal(1,0)*AA(0,0)+SI_tmateqgal(1,1)*AA(1,0)
     &     +SI_tmateqgal(1,2)*AA(2,0)
      BB(1,1)=SI_tmateqgal(1,0)*AA(0,1)+SI_tmateqgal(1,1)*AA(1,1)
     &     +SI_tmateqgal(1,2)*AA(2,1)
      BB(1,2)=SI_tmateqgal(1,0)*AA(0,2)+SI_tmateqgal(1,1)*AA(1,2)
     &     +SI_tmateqgal(1,2)*AA(2,2)
      BB(2,0)=SI_tmateqgal(2,0)*AA(0,0)+SI_tmateqgal(2,1)*AA(1,0)
     &     +SI_tmateqgal(2,2)*AA(2,0)
      BB(2,1)=SI_tmateqgal(2,0)*AA(0,1)+SI_tmateqgal(2,1)*AA(1,1)
     &     +SI_tmateqgal(2,2)*AA(2,1)
      BB(2,2)=SI_tmateqgal(2,0)*AA(0,2)+SI_tmateqgal(2,1)*AA(1,2)
     &     +SI_tmateqgal(2,2)*AA(2,2)
      CC(0,0)=BB(0,0)**2
      CC(0,1)=BB(0,1)**2
      CC(0,2)=BB(0,2)**2
      CC(1,0)=BB(1,0)**2
      CC(1,1)=BB(1,1)**2
      CC(1,2)=BB(1,2)**2
      CC(2,0)=BB(2,0)**2
      CC(2,1)=BB(2,1)**2
      CC(2,2)=BB(2,2)**2
c *** Check Johnson & Soderblom 1987, column or row major? ***
c *** I believe this is correct following the example stars
c *** from the above paper. *** 
      s1=sfac*BB(1,0)*BB(2,0)
      s2=sfac*BB(1,1)*BB(2,1)
      s3=sfac*BB(1,2)*BB(2,2)
c *** End check, effects above three lines ***
      Ue=s1+CC(0,0)*r1+CC(0,1)*r2+CC(0,2)*r3
      Ve=s1+CC(1,0)*r1+CC(1,1)*r2+CC(1,2)*r3
      We=s1+CC(2,0)*r1+CC(2,1)*r2+CC(2,2)*r3
      Ue=dsqrt(Ue)
      Ve=dsqrt(Ve)
      We=dsqrt(We)
c *** vl & vb ***
      c1=dcos(ao(5)*1000.0d0*mas)*dcos(ao(4)*1000.0d0*mas)
      c2=dcos(ao(5)*1000.0d0*mas)*dsin(ao(4)*1000.0d0*mas)
      c3=dsin(ao(5)*1000.0d0*mas)
      l1=SI_tmatgaleq(0,0)*c1+SI_tmatgaleq(1,0)*c2+
     &     SI_tmatgaleq(2,0)*c3
      l2=SI_tmatgaleq(0,1)*c1+SI_tmatgaleq(1,1)*c2+
     &     SI_tmatgaleq(2,1)*c3
      l3=SI_tmatgaleq(0,2)*c1+SI_tmatgaleq(1,2)*c2+
     &     SI_tmatgaleq(2,2)*c3
      Tao4=ao(4)*1000.0d0*mas
      Tao5=ao(5)*1000.0d0*mas
      call equtogal(ao(1),ao(2),lo_p,bo_p)
      if (lo_p.lt.0.0d0) then
         lo_p=lo_p+2.0d0*pi
      endif
      sinphi=0.40893064d0*dsin(pi/2.0d0+ao(1))/(pi/2.0d0-bo_p)
      phi=dasin(sinphi)
c *** Quadrant loop for phi ***
      if (dsin(sinphi).lt.0.0d0) then
         if (dcos(sinphi).gt.0.0d0) then
            phi=2.0d0*pi+dasin(sinphi)
         endif
      endif
      if (dcos(sinphi).lt.0.0d0) then
         phi=pi-dasin(sinphi)
      endif
c *** Convert back to as/yr ***
      vlo_p=(Tao4*dcos(phi)+Tao5*dsin(phi))/1000.0d0/mas
      vbo_p=(Tao5*dcos(phi)-Tao4*dsin(phi))/1000.0d0/mas

      return
      end

      subroutine equtogal(alfa,delta,gl,gb)
      implicit real*8(a-h,o-z)
      real*8 lq,lpl
      include 'const_math.h'
      include 'const_ast.h'
c      INCLUDE 'const_pot.h'

      T6=60.d0
      T36=3600.d0
      RD=PI/180.d0
      LQ=122.92851d0*RD
      ALFAG=(12.d0+51.d0/T6+26.2755d0/T36)*15.d0*RD
      DECG=(27.d0+7.d0/T6+41.704d0/T36)*RD
      SDG=dSIN(DECG)
      CDG=dCOS(DECG)
 
      sd=dsin(delta)
      cd=dcos(delta)

      sb=sdg*sd+cdg*cd*dcos(alfa-alfag)
      cb=sqrt(1.d0-sb**2.d0)
      gb=datan2(sb,cb)
      slpl=(cd*dsin(alfa-alfag))/cb
      clpl=(cdg*sd-sdg*cd*dcos(alfa-alfag))/cb
      lpl=datan2(slpl,clpl)
      if(lpl.le.0.d0)lpl=lpl+2.d0*pi
      gl=lq-lpl

      return
      end
