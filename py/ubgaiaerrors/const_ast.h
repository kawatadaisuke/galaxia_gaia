c----------------------------------------------------------------------
c  This is 'const_ast.h'
c  Astronomical and physical constants (SI units)
c  ep0 is the reference epoch for the astrometric parameters,
c  reckoned from J2000.0
  
      double precision pc,Gte,time
c en cm
      PARAMETER (pc = 3.086d18,
c en  dines cm2/g2 (cm3/(s2 g)) (Bowers)
     &           Gte = 6.670d-8,
c canvi de km/s/kpc a 1/(100 Myr)  
     &           time=0.10225d0)

      double precision kt,a1,a2,Ag,Dg
      parameter (kt=4.7404705d0)
      parameter (a1=(62.87124882d0*deg))
      parameter (a2=(32.93680516d0*deg))
      parameter (Ag=(282.8594813d0*deg))
      parameter (Dg=(27.12825111d0*deg))



