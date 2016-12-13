!+
MODULE Methods 
! -----------------------------------------------------------------------
! PURPOSE - Compute pressure coefficient using various methods of   
!    hypersonic aerodynamics.                                       
! AUTHOR  - Ralph L. Carmichael, Public Domain Aeronautical Software
!  based on original work of Arvel E. Gentry, Douglas N. Smyth and
!  Wayne R. Oliver of Douglas Aircraft.
!  Original work was sponsored by United States Air Force Wright Laboratories
! REVISION HISTORY                                                  
!   DATE  VERS PERSON  STATEMENT OF CHANGES                         
!   1964?   -    AEG   Original program and inspiration             
!   1973    IV   AEG,DNS,WRO  Release of Mark IV version      
! 30Mar88  0.1   RLC   Original coding (in Pascal)                  
! 12May96  0.4   RLC   Converted to a module for Fortran (90)
! 30Jun96  0.5   RLC   All routines are functions, not subroutines
! 07Dec02  0.6   RLC   Brought up-to-date with an unofficial Mark V
! 21Dec02  0.7   RLC   Added NewtonianPrandtlMeyer and VanDykeUnified
! 24Jun07  0.8   RLC   Added TangentConeEdwards

! NOTES - If you ever see fs(6), that is Mach number
!  Throughout the routines, delta is in degrees and deltar is in radians

  IMPLICIT NONE
  
! ---------------------------------------------------------------------------
  
  CHARACTER(LEN=*),PARAMETER:: METHODS_VERSION = "0.8 (24 Dec 2007)"

  REAL,PARAMETER,PRIVATE:: ONE=1.0
  REAL,PARAMETER,PRIVATE:: PI=3.14159265
  REAL,PARAMETER,PRIVATE:: GAMMA = 1.4
  REAL,PARAMETER,PRIVATE:: GAMMA1 = GAMMA+ONE
  REAL,PARAMETER,PRIVATE:: GAMMA12 = 0.5*GAMMA1
  REAL,PARAMETER,PRIVATE:: GAMMAM1 = GAMMA-ONE
  REAL,PARAMETER,PRIVATE:: GAMMAM1R = 1.0/GAMMAM1
  
  REAL,PARAMETER,PRIVATE:: SQRT6 = 2.4494897
  REAL,PARAMETER,PRIVATE:: TERM2 = 2/GAMMA1
  REAL,PARAMETER,PRIVATE:: TERM3 = GAMMAM1/2
  REAL,PARAMETER,PRIVATE:: NUMAX = 2.27685316        ! 0.5*PI*(Sqrt(6)-1) 
  REAL,PARAMETER,PRIVATE:: CON7  = 2.0/GAMMA
  REAL,PARAMETER,PRIVATE:: EXPT  = GAMMA/GAMMAM1
                   
! Constants for Inverse Prandtl-Meyer Calculations
!   from I.M. Hall, Aeronautical Journal, Sept 1975, p.417
  REAL,PARAMETER,PRIVATE:: IPM1 =  1.3604
  REAL,PARAMETER,PRIVATE:: IPM2 =  0.0962
  REAL,PARAMETER,PRIVATE:: IPM3 = -0.5127
  REAL,PARAMETER,PRIVATE:: IPM4 = -0.6722
  REAL,PARAMETER,PRIVATE:: IPM5 = -0.3278

CONTAINS

!+
FUNCTION ACMempirical(mach,deltar) RESULT(cp)
! ---------------------------------------------------------------------------
! PURPOSE - Calculate pressure based on ACM empirical data
  REAL,INTENT(IN):: mach,deltar
  REAL:: cp
  REAL:: delta
  REAL:: msq   ! mach**2
  REAL:: pe
!----------------------------------------------------------------------------
  msq=mach*mach
  delta=(180/PI)*deltar
  cp=delta/(16.0*msq)
  pe=-1.0/msq
  cp=MAX(cp,pe)
  RETURN
END Function ACMempirical   ! -----------------------------------------------  

!+
FUNCTION BlastWave(mach,cpstag,etac,enpm,xcent) RESULT(cp)
! ---------------------------------------------------------------------------
! PURPOSE - Calculate pressure using blast wave analysis
! REFERENCE: J. Lukasiewicz, Hypersonic Flow-Blast Analogy. AEDC-TR-61-4
! NOTES - This option is difficult to use, because you must add this value to
!  the value obtained by an alternate scheme.
! Drop it???
  REAL,INTENT(IN):: mach,cpstag
  REAL,INTENT(IN):: etac,enpm,xcent
  REAL:: cp
  REAL:: msq   ! mach**2
!----------------------------------------------------------------------------
  msq=mach*mach
  IF (cpstag > 0.5) THEN
    cp=(0.121*msq*etac/(enpm-xcent)**.667 + 0.56)/(GAMMA/2.0*msq)
  ELSE
    cp=(0.067*msq*etac/(enpm-xcent) + 0.44 )/ (GAMMA/2.0*msq)
  END IF  
  RETURN
END Function BlastWave   ! --------------------------------------------------  

!+
FUNCTION BluntBodyViscous() RESULT(cp)
! ---------------------------------------------------------------------------
  REAL:: cp
!----------------------------------------------------------------------------
!  CALCULATE BLUNT BODY VISCOUS EFFECTS    [case 8]
!  190 cp = 0.0
!  ivisin = 1
  
!  THE VISCOUS FORCE COEFFICIENT TAU IS CALCULATED IN
!  SUBROUTINE BLUNT, WHICH ONLY NEEDS TO BE CALLED ONCE
!  FOR EACH SECTION.
  
!  IF(l == 1)CALL blunt(pfs,fs(6),tfs,vis,rhofs,etac,reno,tau,ivisin)
!  shear = tau*COS(deltar)
!  GO TO 220
!220 IF (shear > 1.0E-25) GO TO 230
!  shearx = 0.0
!  sheary = 0.0
!  shearz = 0.0
!  GO TO 410
!  230 shearx = shear * sx
!  sheary = shear * sy
!  shearz = shear * sz
!  GO TO 410

  cp=0
  RETURN
END Function BluntBodyViscous   ! -------------------------------------------  

!+
FUNCTION ConeAtAngleOfAttack() RESULT(cp)
! ---------------------------------------------------------------------------
  REAL:: cp
!----------------------------------------------------------------------------
!  CALCULATE PRESSURE USING CONE AT ANGLE OF ATTACK   [case 6]
!  170 CONTINUE
!  JONES METHOD IS USED
!  it = 0
!  CALCULATE LOCATION OF WINDWARD PLANE
!  175 phiw = 0.0
!  IF (vy == 0.0  .AND.  vz >= 0.0)  GO TO 176
!  phiw = 3.1415926536
!  IF (vy == 0.0)  GO TO 176
!  phiw = ATAN2(-vy,vz)
!  176 CONTINUE
  
!  ANGLE OF ATTACK IN WINDWARD PLANE
!  alfwd =ACOS(-vx/fs(7))
  
!  CONE ANGLE
!  dcr = ASIN(nx)
  
!  MERIDIAN ANGLE
!  IF (ny /= 0.0)  GO TO 177
!  phi = 0.0
!  IF (nz == 0.0)  GO TO 178
!  177 phi = ATAN2(ny, -nz)
  
!  LOCATION FROM WINDWARD PLANE
!  178 phit = phi - phiw
  
!  MACH NUMBER IN PHIT PLANE
!  vc = vz*COS(phi)  - vy*SIN(phi)
!  amp = fs(6)*SQRT(vx**2 + vc**2)/fs(7)
  
!  ANGLE OF AMP TO AXIS
!  IF (vc /= 0.0)  GO TO 179
!  alfp = 0.0
!  IF (vx == 0.0)  GO TO 181
!  179 alfp = ATAN2( vc, -vx)
  
  
!  181 angle(1) = dcr/rc
!  IF (iprck == 1) WRITE(tapeot,1000) phiw, phi,phit, alfwd,alfp,dcr,amp
!  1000 FORMAT(1H0, 7X,4HPHIW,11X,3HPHI, 12X,4HPHIT, 11X,5HALFWD,  &
!      10X,4HALFP, 11X,3HDCR, 12X,3HAMP/1H , 7F15.6)
!  CALL acone(angle,cp,alfwd,phit,alfp,amp,it,iprck)
  
!  GO TO 410


  cp=0
  RETURN
END Function ConeAtAngleOfAttack   ! ----------------------------------------

!+
FUNCTION DahlemBuck(mach,deltar,cosdel) RESULT(cp)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the pressure coeff. using the modified Dahlem-Buck method
! ??? different form for expansion ???
  REAL,INTENT(IN):: mach,deltar,cosdel
  REAL:: cp,cpOriginal
  REAL:: a1,a2
  REAL,PARAMETER:: A225 = 22.5*PI/180.0   ! 22.5 degrees in radians
  REAL:: factor
  REAL:: xlnm
!----------------------------------------------------------------------------
  IF (deltar <= 0.0) THEN
    cp=0.0
    RETURN
  END IF
  
  IF (deltar > A225) THEN   ! first compute the original
    cpOriginal=2.0          ! (cp=cpOriginal*cosdel**2 )
  ELSE    
    cpOriginal=1.0/(ABS(SIN(4.0*deltar)))**0.75+1.0
    IF (cpOriginal > 5.0) cpOriginal=5.0
    IF (cpOriginal < 2.0) cpOriginal=2.0
  END IF  
  
! compute low Mach correction factor. Thus makes it the MODIFIED Dahlem-Buck
  xlnm=LOG(mach)
  a1 = (6.0-0.3*mach)+SIN(PI*(xlnm-0.588)/1.20)
  a2 = -1.15 - 0.5*SIN(PI*(xlnm-0.916)/3.29)
  factor=1.0 + a1*(deltar*180/PI)**a2
  cp=cpOriginal*cosdel*cosdel*factor
  RETURN
END Function DahlemBuck   ! -------------------------------------------------

!+
FUNCTION DeltaWingCorrelation(mach, deltar) RESULT(cp)
! ---------------------------------------------------------------------------
! PURPOSE - Calculate pressure using 
  REAL,INTENT(IN):: mach,deltar
  REAL:: cp
  REAL:: deldlw
  REAL:: emns
!----------------------------------------------------------------------------
  deldlw=MAX(0.1745, deltar)  ! at least one degree
  emns = mach* SIN(deldlw)
  emns = 1.09*emns + EXP(-0.49 *emns)
  cp = 1.66667*(emns*emns-1.0)/mach**2
  RETURN
END Function DeltaWingCorrelation   ! ---------------------------------------

!+
FUNCTION HankeyFlatSurface(mach, deltar, cosdel) RESULT(cp)
! ---------------------------------------------------------------------------
! PURPOSE - Calculate pressure using Hankey flat surface empirical correlation 
! Ref: W.L. Hankey: Optimization of Lifting Re-Entry Vehicles. ASD-TDR-62-1102
!  (March 1963)
! NOTE - The equations shown here are from the Mark4 source code, but they do
!  not agree with the manual. Hm??
  REAL,INTENT(IN):: mach,deltar,cosdel
  REAL:: cp
  REAL:: delta   ! deltar in degrees
  REAL:: hankey
!----------------------------------------------------------------------------
  delta=(180/PI)*deltar
  IF (delta < 10.0) THEN
    hankey= (0.195+0.222594/mach**0.3-0.4)*delta + 4.0
  ELSE
    hankey= 1.95 + 0.3925/(mach**0.3*TAN(deltar))
  END IF  

  cp = hankey* cosdel*cosdel
  RETURN
END Function HankeyFlatSurface   ! ------------------------------------------

!+
FUNCTION InversePrandtlMeyer(nu) RESULT(mach)
! ---------------------------------------------------------------------------
! PURPOSE - Inverse Prandtl-Meyer Function. A simple rational polynomial
!  curve fit, good to 5 or 6 significant figures. Refer to the function
!  InversePrandtlMeyerPrecise if you need full double precision accuracy.
! REF: I.M. Hall
  REAL,INTENT(IN):: nu
  REAL:: mach
  REAL:: y
!----------------------------------------------------------------------------
  y=(nu/NUMAX)**0.6666667
  mach=(1.0+y*(IPM1+y*(IPM2+y*IPM3)))/(1.0+y*(IPM4+y*IPM5))
  RETURN
END Function InversePrandtlMeyer   ! ----------------------------------------

FUNCTION InversePrandtlMeyerPrecise(nu) RESULT(mach)
! ---------------------------------------------------------------------------
! PURPOSE - Inverse Prandtl-Meyer function with high precision. Use Hall's
!  approximation for a good first guess, then apply Newton's method to get
!  greater accuracy. Instead of putting a test for convergence in the
!  algorithm, I studied the function for various values of nu and found that
!  four steps will give full convergence to double precision (64 bits). One
!  step would give adequate precision for single precision.
!  Note the use of beta instead of Mach as the dependant variable until the
!  very last step.
  REAL,INTENT(IN):: nu
  REAL:: mach

  INTEGER,PARAMETER:: MAX = 4  !   1 is enough for single; 4 for double
  INTEGER:: i 
  REAL:: beta,betasq 
  REAL:: err 
  REAL,PARAMETER:: ONE = 1.0, SIX = 6.0
!----------------------------------------------------------------------------
  beta=SQRT(ABS((InversePrandtlMeyer(nu))**2 - ONE))            ! first guess
  DO i=1,MAX                                           ! use Newton MAX times
    err=SQRT6*ATAN(beta/SQRT6)-ATAN(beta)-nu                    ! error in nu 
    betasq=beta*beta
    beta=beta-err*(SIX+betasq)*(ONE+betasq)/(SIX*TERM2*betasq)  ! d(nu)/d(beta)
  END DO
  mach=SQRT(beta*beta+ONE)
  RETURN
END Function InversePrandtlMeyerPrecise   ! ---------------------------------

!+
FUNCTION Newtonian(cpstag,deltar) RESULT(cp)
! ---------------------------------------------------------------------------
! PURPOSE - This is really modified Newtonian because of the option of 
!  choosing something other than 2.0 as stagnation pressure coefficient
  REAL,INTENT(IN):: cpstag,deltar
  REAL:: cp
!----------------------------------------------------------------------------
  cp=cpStag*Sin(deltar)**2
  RETURN
END Function Newtonian   ! --------------------------------------------------

!+
FUNCTION NewtonianPrandtlMeyer(mach,deltar,cpstag) RESULT(cp)
! ---------------------------------------------------------------------------
! PURPOSE - 

! REFERENCE - L.G. Kaufman: Pressure Estimation Techniques for Hypersonic
!  Flows over Blunt Bodies. Journal of the Astronautical Sciences.
!  Vol. X, No. 2, Summer 1963.
  REAL,INTENT(IN):: mach,deltar,cpstag
  REAL:: cp
  REAL:: emlow,emup   ! limits for matching point Mach
  INTEGER:: k
  REAL:: m1,m2
  INTEGER,PARAMETER:: MAX_ITER = 20
  REAL:: mdelta
  REAL:: msq
  REAL:: msubq
  REAL:: nu
  REAL:: p1,p2
  REAL:: pc
  REAL:: pcap  ! ratio of freestream static to freestream stagnation
  REAL:: ppfs
  REAL:: ppo
  
  REAL:: q  ! ratio of matching point pressure to freestream static
  REAL:: sdeltq
!----------------------------------------------------------------------------
  msq=mach*mach
  pcap=(TERM2/msq)**EXPT * ((2.0*GAMMA*msq-GAMMA+1.0)/GAMMA1)**GAMMAM1R
  
!... calculate matching Mach number by iteration
  m2=0.0
  p2=0.0
  emlow=0.91+0.3125*GAMMA
  emup=emlow+0.4
  msubq=emlow
  DO k=1,MAX_ITER
    q=(2.0/(2.0+GAMMAM1*msubq**2))**EXPT
    pc=q*(1.0-(GAMMA**2*msubq**4*q)/(4.0*(1.0-q)*(msubq**2-1.0)))
    IF (ABS(msubq-m2) < 1E-4) EXIT   ! success
    p1=p2
    p2=pc
    m1=m2
    m2=msubq
    IF (k==1) THEN
      msubq=emup
      CYCLE
    END IF
    IF (ABS(p2-p1) < 1E-6) EXIT   ! success
    msubq=m1+(pcap-p1)*(m2-m1)/(p2-p1)
    msubq=MIN(msubq,emup)   ! don't let msubq get outside the limits
    msubq=MAX(msubq,emlow)  
  END DO

!... if flow has not reached the matching point, forget it and use Newtonian
  IF (q < pcap) THEN
    cp=cpStag*SIN(deltar)**2
    RETURN
  END IF  

!... calculate surface slope at matching point
  sdeltq=SQRT((q-pcap)/(1.0-pcap))
  
! get Mach number from the inverse PrandtlMeyer function
  nu=ASIN(sdeltq)+deltar
  mdelta=InversePrandtlMeyer(nu)

! calculate the surface pressure ratio
  ppo=(1.0+0.5*GAMMAM1*mdelta**2)**(-EXPT)

! calculate the surface to freestream pressure ratio
  ppfs=ppo/pcap

! calculate the surface pressure coefficient
  cp=(2.0/(GAMMA*mach**2))*(ppfs-1.0)
  RETURN
END Function NewtonianPrandtlMeyer   ! -------------------------------------- 

!+
FUNCTION OSUBluntBody(cpStag,mach,deltar) RESULT(cp)
! ---------------------------------------------------------------------------
! PURPOSE - 
! NOTE - This seems to be a left-over method from previous versions of the
!  program. Reference unknown. Use at your own risk.
  REAL,INTENT(IN):: cpStag,mach,deltar
  REAL:: cp

  REAL:: theta, xx, xy 
!----------------------------------------------------------------------------
  theta=PI/2-deltar
  xx=0.32+0.455*Cos(theta)     + 0.195*Cos(2.0*theta)+     &
          0.035*Cos(3.0*theta) - 0.005*Cos(4.0*theta)
  xy=0.7*mach**2
  cp=(xx*(cpStag*xy+1.0)-1.0)/xy
  RETURN
END Function OSUBluntBody   ! -----------------------------------------------

!+
FUNCTION PrandtlMeyerXXX(mach, deltar) RESULT(cp)
! ---------------------------------------------------------------------------
  REAL,INTENT(IN):: mach,deltar
  REAL:: cp

  INTEGER,PARAMETER:: MAX = 5
  REAL:: machsq, betasq, beta, nu, del, bracket 
  REAL:: nuZero,err 
  INTEGER:: i 
!----------------------------------------------------------------------------
  machsq=mach*mach
  beta=Sqrt(Abs(machsq-1.0))
  nuZero=SQRT6*ATan(beta/SQRT6)-ATan(beta) !   free-stream nu }
  del=-deltar !  should be > 0 }
  nu=nuZero+del !      nu on the expansion surface }
  beta=5.0/(NUMAX-nu) !   make a guess for beta - infinite Mach approx. }
  beta=Sqrt(Sqrt(nu/NUMAX))*5.0/(NUMAX-nu)
  err=SQRT6*ATan(beta/SQRT6)-ATan(beta)-nu !   error in nu at this beta}
  DO i=1,MAX  ! use Newton-Raphson MAX times 
    betasq=beta**2
    beta=beta-err*(1+betasq/6.0)*(1+betasq)/(TERM2*betasq)       !  estimate 
    err=SQRT6*ATAN(beta/SQRT6)-ATAN(beta)-nu !   new error in nu }
  END DO
  bracket=(1.0+TERM3*(1+beta*beta))/(1.0+TERM3*machsq)
  cp=(CON7/machsq)*(bracket**(-EXPT) -1.0)
  RETURN
END Function PrandtlMeyerXXX   ! --------------------------------------------

!+
FUNCTION PrandtlMeyer(mach,deltar) RESULT(cp)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the pressure coefficient associated with an change of
!  freestream direction of deltar. Since deltar is negative in the calling
!  program, note that you subtract deltar from nuZero to get the nu on the
!  expansion surface.
  REAL,INTENT(IN):: mach
  REAL,INTENT(IN):: deltar  ! radians  ( will be negative)
  REAL:: cp

  REAL:: machsq, localMachsq, beta, nu, bracket,nuZero
!----------------------------------------------------------------------------
  machsq=mach*mach
  beta=Sqrt(Abs(machsq-1.0))
  nuZero=SQRT6*ATan(beta/SQRT6)-ATan(beta)                   ! free-stream nu
  nu=nuZero-deltar                              ! nu on the expansion surface
  IF (nu > NUMAX) THEN
    cp=-2.0/(GAMMA*machsq)   ! can't go any lower than this
    RETURN
  END IF  
  localMachsq=InversePrandtlMeyer(nu)**2
  bracket=(1.0+TERM3*localMachsq)/(1.0+TERM3*machsq)
  cp=(CON7/machsq)*(bracket**(-EXPT)-1.0)
  RETURN
END Function PrandtlMeyer   ! -----------------------------------------------

!+
FUNCTION SmythDeltaWing(mach,deltar) RESULT(cp)
! ---------------------------------------------------------------------------
  REAL,INTENT(IN):: mach,deltar
  REAL:: cp
  REAL:: deldlw
  REAL:: emns
!----------------------------------------------------------------------------
  deldlw = MAX(0.01745, deltar)      ! nothing less than one degree
  emns=mach*SIN(deldlw)
  emns=1.09*emns + EXP(-0.49*emns)   ! empirical 
  cp=1.66667*(emns**2-1.0)/mach**2
  RETURN
END Function SmythDeltaWing   ! ---------------------------------------------

!+
FUNCTION OldTangentCone(mach,deltar) RESULT(cp)
! ---------------------------------------------------------------------------
! This is the old tangent-cone empirical method, described on pp. 125-128 of
!  volume 2 of AFFDL-TR-73-159. This was the primary method in the Mark3
!  program but is superceded in Mark4 and Mark5.
! It is not as accurate, but it is simple and fast.
  REAL,INTENT(IN):: mach,deltar
  REAL:: cp

  REAL:: s,xx 
!----------------------------------------------------------------------------
  s=Sin(deltar)
  xx=mach*s
  xx=1.090909*xx+Exp(-0.5454545*xx)
  xx=xx*xx
  cp=48.0*xx*s*s/(23.0*xx-5.0)
  RETURN
END Function OldTangentCone   ! ---------------------------------------------

!+
FUNCTION TangentConeEdwards(mach,deltar) RESULT(cp)
! ------------------------------------------------------------------------------
! Ref: NASA TP 1539 (Appendix)
  REAL,INTENT(IN):: mach, deltar
  REAL:: cp
  REAL:: machNs,m2   ! mach number normal to the shock and its square
  REAL:: s   ! sin(deltar)
!-------------------------------------------------------------------------------
  s=SIN(deltar)
  machNs=(0.87*mach-0.544)*s + 0.53
  m2=machNs**2
  cp=48.0*m2*s**2/(23.0*m2-5.0)
  RETURN
END Function TangentConeEdwards   ! --------------------------------------------

!+
FUNCTION TangentWedgeEmpirical(mach, deltar) RESULT(cp)
! ---------------------------------------------------------------------------
  REAL,INTENT(IN):: mach, deltar
  REAL:: cp

  REAL:: xx 
!----------------------------------------------------------------------------
  xx=mach*SIN(deltar)
  cp=((1.2*xx+EXP(-0.6*xx))**2-1.0)/(0.6*mach**2)
  RETURN
END Function TangentWedgeEmpirical   ! --------------------------------------

!+
FUNCTION TangentWedge(mach, deltar) RESULT(cp)
! ---------------------------------------------------------------------------
  REAL,INTENT(IN):: mach,deltar
  REAL:: cp

  REAL:: machsq, mach4, sn2, b, c, d
  REAL:: Q,R, disc, costh, theta 
  REAL:: root1, root2, root3 
!----------------------------------------------------------------------------
! If deltar is greater than 45.585 degrees, then the shock is detached, 
!  regardless of the value of mach. Use TangentWedgeEmpirical
  IF (deltar > 0.7956) THEN
    cp=TangentWedgeEmpirical(mach,deltar)
    RETURN
  END IF
    
  machsq=mach**2
  mach4=machsq**2
  sn2=SIN(deltar)**2
  
!... There can be numerical problems with very small wedge angles. Use the
!    equation on p. 92 of Liepmann and Roshko, Gasdynamics.  
  IF (deltar < 0.035) THEN
    cp=GAMMA*machsq*deltar/SQRT(machsq**2-1.0)
  END IF
  
  b=-(machsq+2)/machsq-GAMMA*sn2
  c=(2*machsq+1)/mach4+sn2*(GAMMA12**2+GAMMAM1/machsq)
  d=(sn2-1)/mach4
  q=(b*b-3*c)/9.0
  IF (q==0.0) WRITE(*,*) "Q=0"
  R=(b*(2.0*b**2-9.0*c)+27.0*d)/54.0
  disc=q*q*q - r*r
  IF (disc < 0.0) THEN
    cp=TangentWedgeEmpirical(mach,deltar)
  ELSE
    costh=r/SQRT(q*q*q)
    theta=ACOS(costh)
    c=-2*SQRT(q)
    d=b/3
    root1=c*COS(theta/3)-d
    root2=c*COS((theta+2*PI)/3)-d
    root3=c*COS((theta+4*PI)/3)-d
    ! root3 is supposed to be the middle root, but check it anyway
    IF (root1>root3) WRITE(*,*) 'VIOLATION #1 ', mach, deltar, root1,root2,root3
    IF (root2<root3) WRITE(*,*) 'VIOLATION #2 ', mach, deltar, root1,root2,root3
    IF (root3 < 0.0) WRITE(*,*) "VIOLATION #3 ", mach, deltar, root1,root2,root3
    IF (ABS(root3) > 1.0) WRITE(*,*) "VIOLATION #4", mach, deltar, root1,root2,root3
    theta=ASIN(SQRT(root3))   ! wave angle 
    cp=4.0*(machsq*root3-1.0)/(machsq*(GAMMA+1.0))
  END IF
  RETURN
END Function TangentWedge   ! -----------------------------------------------

!+
FUNCTION TangentWedgeInfiniteMach(mach,deltar) RESULT(cp)
! ---------------------------------------------------------------------------
! PURPOSE - Calculate pressure using tangent wedge - infinite Mach method
!    add theory


  REAL,INTENT(IN):: mach,deltar
  REAL:: cp

  REAL:: emns   ! ???
!----------------------------------------------------------------------------
  emns = 0.5*GAMMA1*mach*SIN(deltar)+EXP(-0.25*GAMMA1*mach*SIN(deltar))
  cp = (4./GAMMA1)*(emns**2-1.)/mach**2
  RETURN
END Function TangentWedgeInfiniteMach   ! -----------------------------------  

!+
FUNCTION VanDykeUnified(mach, deltar) RESULT(cp)
! ---------------------------------------------------------------------------
! PURPOSE - Calculate pressure coefficient using the Van Dyke unified theory.
! REFERENCE: Milton Van Dyke: A Study of Hypersonic Small Disturbance Theory.
!  NACA Report 1194, 1954.
! NOTE - This does not agree !!!!
  REAL,INTENT(IN):: mach,deltar
  REAL:: betasq
  REAL:: bracket
  REAL:: cp
  REAL:: cpVac
  REAL:: gammaTerm
  REAL::h,hsq  ! hypersonic similarity parameter (and its square)
  REAL:: machsq   ! square of the Mach number
!----------------------------------------------------------------------------
  IF (deltar==0.0) THEN
    cp=0.0
    RETURN
  END IF
    
  machsq=mach*mach
  betasq=machsq-1.0
  hsq=betasq*deltar**2
  h=SQRT(hsq)
  gammaTerm=0.5*GAMMA1  
  
  IF (deltar > 0.0) THEN
    cp=deltar**2*(gammaTerm + SQRT(gammaTerm**2+4.0/hsq))
  ELSE
    bracket=(1.0-0.5*GAMMAM1*h)**(2.0*EXPT) - 1.0
    cp=2.0*deltar**2*bracket/(GAMMA*hsq) 
    cpVac= -2.0/(GAMMA*machsq) 
    cp=MAX(cp,cpVac)
  END IF

  RETURN
END Function VanDykeUnified   ! ---------------------------------------------


END Module Methods   ! ======================================================
