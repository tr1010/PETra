!+
INCLUDE 'methods.f90'
!+
MODULE NetworkProcedures
! ---------------------------------------------------------------------------
! PURPOSE - Define the derived type Network and collect the procedures that
!  operate on data with the type Network. 

! AUTHORS  - Ralph L. Carmichael, Public Domain Aeronautical Software
!            Charlotte Craidon, NASA Langley, coordinator of the Langley
!               wireframe geometry standard as documented in NASA TM 85767. 

! REVISION HISTORY                                                  
!   DATE  VERS PERSON  STATEMENT OF CHANGES   
!   1977         CC    Published NASA TM 85767
!   1989   0.1   RLC   Various routines to read networks in LaWgs format
! 23Nov02  0.5   RLC   Rewrote completely (Fortran 95)
! 07Dec02  0.6   RLC   First working version of the hypersonic aerodynamics program


IMPLICIT NONE

  TYPE:: Network
    CHARACTER(LEN=80):: name
    INTEGER:: id
    INTEGER:: compMethod, expMethod
    INTEGER:: rows,cols
    LOGICAL:: reversed
    INTEGER:: symLocal, symGlobal
    INTEGER:: localImage, globalImage ! index of original network; =0 if an original
    REAL,POINTER,DIMENSION(:,:):: x,y,z
    REAL:: fx,fy,fz, mx,my,mz   ! total force and moment on the network
    REAL,DIMENSION(3):: rotate,translate,scale
    REAL,DIMENSION(3):: averageCenter, averageNormal   ! averages
    REAL:: totalArea
  END TYPE Network
  

  INTEGER,PARAMETER,PRIVATE:: DBG = 3
  CHARACTER(LEN=*),PARAMETER:: NETWORK_VERSION = '0.64 (27Dec02)'
  
  PRIVATE:: BuildRotationMatrix
  PUBLIC::  CreateGlobalImageNetworks
  PUBLIC::  CreateLocalImageNetworks  
  PUBLIC::  DeallocateNetworks
  PUBLIC::  DumpNetworkGridPoints
  PUBLIC::  DumpNetworkProperties
  PUBLIC::  ReadOneNetwork
  PUBLIC::  ScanNetworks
  PRIVATE:: StripFortranQuotes
! ---------------------------------------------------------------------------

CONTAINS


!+
FUNCTION BuildRotationMatrix(r) RESULT(rot)
! ---------------------------------------------------------------------------
! PURPOSE - Form the rotation matrix used for geometrical transformations
!  Taken from NASA TM 85767 defining LaWGS

  REAL,INTENT(IN),DIMENSION(:):: r   ! rotation angles, degrees
  REAL,DIMENSION(3,3):: rot

  REAL,PARAMETER:: PI=3.14159265, RAD=PI/180

  REAL:: cPhi,cTheta,cPsi
  REAL:: sPhi,sTheta,sPsi
!----------------------------------------------------------------------------
  cPhi=COS(RAD*r(1))
  sPhi=SIN(RAD*r(1))
  cTheta=COS(RAD*r(2))
  sTheta=SIN(RAD*r(2))
  cPsi=COS(RAD*r(3))
  sPsi=SIN(RAD*r(3))

  rot(1,1)= cTheta*cPsi
  rot(2,1)= cTheta*sPsi
  rot(3,1)=-sTheta

  rot(1,2)=-sPsi*cPhi + sTheta*cPsi*sPhi
  rot(2,2)= cPsi*cPhi + sTheta*sPsi*sPhi
  rot(3,2)= cTheta*sPhi

  rot(1,3)= sPsi*sPhi + sTheta*cPsi*cPhi
  rot(2,3)=-cPsi*sPhi + sTheta*sPsi*cPhi
  rot(3,3)= cTheta*cPhi

  RETURN
END FUNCTION BuildRotationMatrix   ! ----------------------------------------

!+
SUBROUTINE ConvertNetwork(a, r, t, s) 
!   -------------------------------------------------------------------------
! PURPOSE - Transform the gridpoints in a network according to the
!   transformation equations. From the LaWgs document, NASA TM-85767, the
!   order of conversion is rotation, translation, scale

  TYPE(Network),INTENT(IN OUT):: a
  REAL,INTENT(IN),DIMENSION(:):: r   ! rotation angles (deg)
  REAL,INTENT(IN),DIMENSION(:):: t   ! translation distances
  REAL,INTENT(IN),DIMENSION(:):: s   ! scale factors

  REAL,DIMENSION(3):: before,after
  INTEGER:: i,j
  REAL,DIMENSION(3,3):: rot
!----------------------------------------------------------------------------
  rot=BuildRotationMatrix(r)

  DO j=1,a%cols
    DO i=1,a%rows
      before(1)=a%x(i,j)
      before(2)=a%y(i,j)
      before(3)=a%z(i,j)
      after=MATMUL(rot,before)
      after=s*(after + t)
      a%x(i,j)=after(1)
      a%y(i,j)=after(2)
      a%z(i,j)=after(3)
    END DO  
  END DO
  RETURN
END Subroutine ConvertNetwork   ! -------------------------------------------

!+
SUBROUTINE CreateGlobalImageNetworks(n,a)
! ---------------------------------------------------------------------------
! PURPOSE - If any networks in the original file have the symGlobal flag = 1,
!  make an image network of it. Note that n and a will be modified and that a
!  must be dimensioned large enough to hold the new networks. Fatal if not.

  INTEGER,INTENT(IN OUT):: n
  TYPE(Network),INTENT(IN OUT),DIMENSION(:):: a 
 
  INTEGER:: i,j,k
  INTEGER:: nOriginal
! ---------------------------------------------------------------------------
  nOriginal=n
  DO k=1,nOriginal
    IF (a(k)%symGlobal==0) CYCLE
    IF (n+1 > SIZE(a)) THEN
      WRITE(*,*) 'Too many networks (2)'
      STOP
    END IF
    n=n+1
    a(n)%name='global image of '//a(k)%name
    a(n)%id=1000+a(k)%id
    a(n)%compMethod=a(k)%compMethod
    a(n)%expMethod=a(k)%expMethod
    a(n)%rows=a(k)%rows
    a(n)%cols=a(k)%cols
    a(n)%reversed = .NOT.a(k)%reversed   ! think about that one
    a(n)%symLocal=a(k)%symLocal   ! but it won't be used
    a(n)%symGlobal=1   ! but it won't be used
    a(n)%localImage=0 
    a(n)%globalImage=k  ! keep track of where it came from
    i=a(n)%rows
    j=a(n)%cols
    ALLOCATE(a(n)%x(i,j), a(n)%y(i,j), a(n)%z(i,j) )
    a(n)%x=a(k)%x
    a(n)%y=-a(k)%y   ! this is where we do the imaging
    a(n)%z=a(k)%z
  END DO
  RETURN
END Subroutine CreateGlobalImageNetworks   ! ---------------------------------

!+
SUBROUTINE CreateLocalImageNetworks(n,a)
! ---------------------------------------------------------------------------
! PURPOSE - If any networks in the original file have the symLocal flag = 1,
!  make an image network of it. Note that n and a will be modified and that a
!  must be dimensioned large enough to hold the new networks. Fatal if not.

  INTEGER,INTENT(IN OUT):: n
  TYPE(Network),INTENT(IN OUT),DIMENSION(:):: a 
 
  INTEGER:: i,j,k
  INTEGER:: nOriginal
! ---------------------------------------------------------------------------
  nOriginal=n
  DO k=1,nOriginal
    IF (a(k)%symLocal==0) CYCLE
    IF (n+1 > SIZE(a)) THEN
      WRITE(*,*) 'Too many networks (1)'
      STOP
    END IF
    n=n+1
    a(n)%name='local image of '//a(k)%name
    a(n)%id=1000+a(k)%id
    a(n)%compMethod=a(k)%compMethod
    a(n)%expMethod=a(k)%expMethod
    a(n)%rows=a(k)%rows
    a(n)%cols=a(k)%cols
    a(n)%reversed = .NOT.a(k)%reversed   ! think about that one
    a(n)%symLocal=1   ! but it won't be used
    a(n)%symGlobal=a(k)%symGlobal   ! we still may need to image it again!
    a(n)%localImage=k   ! keep track of where it came from
    a(n)%globalImage=0   ! will be set by CreateGlobalImageNetworks
    i=a(n)%rows
    j=a(n)%cols
    ALLOCATE(a(n)%x(i,j), a(n)%y(i,j), a(n)%z(i,j) )
    a(n)%x=a(k)%x
    a(n)%y=-a(k)%y   ! this is where we do the imaging
    a(n)%z=a(k)%z
    a(n)%rotate=a(k)%rotate
    a(n)%translate=a(k)%translate
    a(n)%scale=a(k)%scale
  END DO
  RETURN
END Subroutine CreateLocalImageNetworks   ! ---------------------------------

!+
SUBROUTINE DeAllocateNetworks(a)
! ---------------------------------------------------------------------------
! PURPOSE - Deallocate the memory allocated for gridpoints in each network
  TYPE(Network),INTENT(IN OUT),DIMENSION(:):: a  
  INTEGER:: k
!----------------------------------------------------------------------------
  DO k=1,SIZE(a)
    DEALLOCATE(a(k)%x, a(k)%y, a(k)%z)
  END DO
  RETURN
END Subroutine DeAllocateNetworks   ! ---------------------------------------

!+
SUBROUTINE DumpNetworkForcesAndMoments(efu, a)
! ---------------------------------------------------------------------------
! PURPOSE - Print the forces and moments on the defined networks
  INTEGER,INTENT(IN):: efu
  TYPE(Network),INTENT(IN),DIMENSION(:):: a
  
  CHARACTER(LEN=*),PARAMETER:: FMT = '(I4,2I3,6F11.4)'
  INTEGER:: n
!----------------------------------------------------------------------------
  WRITE(efu,*) 'NETWORK FORCES AND MOMENTS'
  WRITE(efu,'(T4,A)') &
   '#  C  E       fx         fy         fz         mx         my         mz'
  DO n=1,SIZE(a)
    WRITE(efu,FMT) n, a(n)%compMethod, a(n)%expMethod, &
      a(n)%fx, a(n)%fy, a(n)%fz, a(n)%mx, a(n)%my, a(n)%mz
  END DO        
  RETURN
END Subroutine DumpNetworkForcesAndMoments   ! ------------------------------

!+
SUBROUTINE DumpNetworkGridPoints(efu, a)
! ---------------------------------------------------------------------------
! PURPOSE - Print the defined networks
  INTEGER,INTENT(IN):: efu
  TYPE(Network),INTENT(IN),DIMENSION(:):: a
  
  CHARACTER(LEN=*),PARAMETER:: FMT = '(2I4,3F12.6)'
  INTEGER:: i,j,n
!----------------------------------------------------------------------------
  DO n=1,SIZE(a)
    WRITE(efu,*) 'NETWORK',n
    WRITE(efu,'(T4,A)') 'R   C       x           y             z'
    DO j=1,a(n)%cols
      DO i=1,a(n)%rows
        WRITE(efu,FMT) i,j, a(n)%x(i,j), a(n)%y(i,j), a(n)%z(i,j)
      END DO
    END DO    
    WRITE(efu,*)   ! blank line
  END DO
  RETURN
END Subroutine DumpNetworkGridPoints   ! ------------------------------------

!+
SUBROUTINE DumpNetworkProperties(efu, a)
! ---------------------------------------------------------------------------
! PURPOSE - Print the defined networks
  INTEGER,INTENT(IN):: efu
  TYPE(Network),INTENT(IN),DIMENSION(:):: a
  
  CHARACTER(LEN=*),PARAMETER:: FMT1 = '(I4,F12.5, 3F12.5,3X,A)'
  CHARACTER(LEN=*),PARAMETER:: FMT2 = '(I4,3F9.4,3X,A)'
  INTEGER:: n
!----------------------------------------------------------------------------
  WRITE(efu,*)
  WRITE(efu,*) 'NETWORK AREAS AND CENTERS'
  WRITE(efu,'(T4,A)') '#   area   x-center   y-center   z-center   name'
  DO n=1,SIZE(a)
    WRITE(efu,FMT1) n, &
      a(n)%totalArea, a(n)%averageCenter(:), Trim(a(n)%name)
  END DO      
  
  WRITE(efu,*)
  WRITE(efu,*) 'NETWORK NORMALS'
  WRITE(efu,'(T4,A)') '#    nx       ny       nz      name'
  DO n=1,SIZE(a)
    WRITE(efu,FMT2) n, a(n)%averageNormal(:), Trim(a(n)%name)
  END DO
  
  WRITE(efu,*)
  WRITE(efu,*) 'NETWORK SYMMETRY CODES AND IMAGE NUMBERS'
  WRITE(efu,*) '   SL=local symmetry code;   SG=global symmetry code'
  WRITE(efu,*) '   LI=network number that this net is a local image of'
  WRITE(efu,*) '   GI=network number that this net is a global image of'
  WRITE(efu,'(T4,A)') '#  SL  SG  LI  GI    name'
  DO n=1,SIZE(a)
    WRITE(efu,'(5I4,3X,A)') n,a(n)%symLocal, a(n)%symGlobal, &
      a(n)%localImage, a(n)%globalImage, Trim(a(n)%name)
  END DO
      
  RETURN
END Subroutine DumpNetworkProperties   ! ------------------------------------

!+
SUBROUTINE PrintNetworkMethodsAndNormals(efu, a)
! ---------------------------------------------------------------------------
! PURPOSE - Print the defined networks
  INTEGER,INTENT(IN):: efu
  TYPE(Network),INTENT(IN),DIMENSION(:):: a
  
  CHARACTER(LEN=*),PARAMETER:: FMT1 = '(I4,2I3,3F7.3,3X,A)'
  INTEGER:: n
!----------------------------------------------------------------------------
  WRITE(efu,*) 'COMPUTATIONAL METHODS AND AVERAGE NORMALS'
  WRITE(efu,*) '   C=compression method;   E=expansion method'
  WRITE(efu,'(T4,A)') '#  C  E    nx     ny     nz'
  DO n=1,SIZE(a)
    WRITE(efu,FMT1) n, &
      a(n)%compMethod, a(n)%expMethod, a(n)%averageNormal(:), Trim(a(n)%name)
  END DO   
  WRITE(efu,*)   ! blank line     
  RETURN
END Subroutine PrintNetworkMethodsAndNormals   ! ----------------------------

!+
SUBROUTINE ReadOneNetwork(efu, a)
! ---------------------------------------------------------------------------
! PURPOSE - Read one network in LaWgs format. Allocate the memory for 
!  variables x,y,z and load them with the grid point coordinates.
  INTEGER,INTENT(IN):: efu
  TYPE(NETWORK),INTENT(OUT):: a
  
  CHARACTER(LEN=132):: buffer
  INTEGER:: errCode
  INTEGER:: i,j,k,n
  INTEGER:: nobj
  CHARACTER(LEN=132):: title
  INTEGER:: rows,cols
  INTEGER::symLocal,symGlobal
  REAL,DIMENSION(3):: rotate,translate,scale
  
  REAL,ALLOCATABLE,DIMENSION(:):: inpX,inpY,inpZ 
!----------------------------------------------------------------------------
  READ(efu, '(A)') title
  CALL StripFortranQuotes(title)
  WRITE(*,*) "Reading network "//Trim(title)
  
  rotate=0.0
  translate=0.0
  scale=1.0
  READ(efu,'(A)') buffer
  READ(buffer,*,IOSTAT=errCode) &
    nobj, cols, rows, symLocal, rotate,translate,scale, symGlobal
  IF (errCode /= 0) THEN
    WRITE(DBG,*) 'Some parameters on the col,row record missing'  
    WRITE(DBG,*) 'buffer:'//Trim(buffer)
  END IF  
  n=rows*cols
  ALLOCATE(inpX(n),inpY(n),inpZ(n))
  Read(efu,*) (inpX(k), inpY(k), inpZ(k), k=1,n)

  ALLOCATE(a%x(rows,cols), a%y(rows,cols), a%z(rows,cols) )

  a%name = title
  a%id = nobj
  a%cols = cols
  a%rows = rows
  a%symLocal=symLocal
  a%symGlobal=symGlobal
  a%rotate=rotate
  a%translate=translate
  a%scale=scale
  a%localImage=0
  a%globalImage=0

  k=0
  DO j=1,cols
    DO i=1,rows
      k=k+1
      a%x(i,j)=inpX(k)
      a%y(i,j)=inpY(k)
      a%z(i,j)=inpZ(k)
    END DO
  END DO    
      
  DEALLOCATE(inpX, inpY, inpZ)
  
  RETURN
END Subroutine ReadOneNetwork   ! -------------------------------------------

!+
SUBROUTINE ScanNetworks(efu, nets,ngrid,npanel,maxRows,maxCols,maxGrid)
! ---------------------------------------------------------------------------
! PURPOSE - Read the Wgs file that is already open as unit efu to determine
!  the number of elements defined by the file. Read networks until 
!  end-of-file or a I/O error is encountered or a net with zero rows or cols.
!  You will probably get warning messages that nobj, x,y,z are
!  never used. True, but we have to read past them in the file. 
  INTEGER,INTENT(IN):: efu
  INTEGER,INTENT(OUT),OPTIONAL:: nets,ngrid,npanel
  INTEGER,INTENT(OUT),OPTIONAL:: maxRows,maxCols,maxGrid

  INTEGER:: k 
  INTEGER:: nobj 
  CHARACTER(LEN=132):: title, buffer
  INTEGER:: rows,cols
  INTEGER:: errCode
  REAL:: x,y,z
!----------------------------------------------------------------------------
  IF (Present(nets))    nets=0
  IF (Present(ngrid))   ngrid=0
  IF (Present(npanel))  npanel=0
  IF (Present(maxRows)) maxRows=0
  IF (Present(maxCols)) maxCols=0
  IF (Present(maxGrid)) maxGrid=0

  REWIND(efu)
  READ(efu, '(A)', IOSTAT=errCode) title     ! general title
  IF (errCode /= 0) WRITE(*,*) "Unable to read general title"

  DO
    READ(efu, '(A)', IOSTAT=errCode) title   ! network title
    IF (errCode /= 0) EXIT
    CALL StripFortranQuotes(title)
    WRITE(*,*) "Scanning network "//Trim(title)
    READ(efu, '(A)', IOSTAT=errCode) buffer
    IF (errCode /= 0) EXIT
    READ(buffer, *, IOSTAT=errCode) nobj, cols, rows
    IF (errCode /= 0) EXIT
    IF (rows <= 0) EXIT 
    IF (cols <= 0) EXIT
    READ(efu,*,IOSTAT=errCode) (x,y,z, k=1,rows*cols)

    IF (Present(nets))    nets=nets+1
    IF (Present(ngrid))   ngrid=ngrid+rows*cols
    IF (Present(npanel))  npanel=npanel+(rows-1)*(cols-1) 
    IF (Present(maxRows)) maxRows=MAX(rows,maxRows)
    IF (Present(maxCols)) maxCols=MAX(cols,maxCols)
    IF (Present(maxGrid)) maxGrid=MAX(rows*cols,maxGrid)
  END DO

  REWIND(efu)
  RETURN
END Subroutine ScanNetworks   ! ---------------------------------------------

!+
SUBROUTINE StripFortranQuotes(t)
! ---------------------------------------------------------------------------
! PURPOSE - If the first non-blank character and the last non-blank
!   character in a string are both APOSTROPHE, then change them to
!   blanks and adjust the resulting string to the left.

  CHARACTER(LEN=*),INTENT(IN OUT):: t

  CHARACTER(LEN=1),PARAMETER:: APOSTROPHE = "'"
!!!  CHARACTER(LEN=LEN(t)):: s   ! elegant, works on Lahey, but not Absoft
  CHARACTER(LEN=256):: s
  INTEGER:: k
!----------------------------------------------------------------------------
  s=Trim(AdjustL(t))   ! get rid of leading and trailing blanks
  k=Len_Trim(s)
  IF (k==0) RETURN   ! just forget it
  IF( (s(k:k) == APOSTROPHE) .AND. (s(1:1)==APOSTROPHE) ) THEN
    s(k:k)=" "
    s(1:1)=" "
  END IF
  t = Trim(AdjustL(s))
  RETURN
END Subroutine StripFortranQuotes   ! ---------------------------------------

!+
SUBROUTINE TransformNetworks(a) 
!   -------------------------------------------------------------------------
! PURPOSE - Transform the gridpoints in an array of networks according to
!  the transformation equations. From the LaWgs document, NASA TM 85767,
!  the order of conversion is rotation, translation, scale.

  TYPE(Network),INTENT(IN OUT),DIMENSION(:):: a

  REAL,DIMENSION(3):: before,after
  INTEGER:: i,j,k
  REAL,DIMENSION(3,3):: rot
!----------------------------------------------------------------------------
  DO k=1,SIZE(a)
    IF (MAXVAL(ABS(a(k)%rotate))==0.0    .AND.   &
        MAXVAL(ABS(a(k)%translate))==0.0 .AND.   &
        MAXVAL(ABS(a(k)%scale-1.0))==0.0 ) CYCLE
        
    rot=BuildRotationMatrix(a(k)%rotate)

    DO j=1,a(k)%cols
      DO i=1,a(k)%rows
        before(1)=a(k)%x(i,j)
        before(2)=a(k)%y(i,j)
        before(3)=a(k)%z(i,j)
        after=MATMUL(rot,before)       ! first rotate (about origin)
        after=after + a(k)%translate   ! then translate
        after=after * a(k)%scale       ! then scale
        a(k)%x(i,j)=after(1)
        a(k)%y(i,j)=after(2)
        a(k)%z(i,j)=after(3)
      END DO  
    END DO
    
  END DO  
  RETURN
END Subroutine TransformNetworks   ! ----------------------------------------

END Module NetworkProcedures   ! ============================================



!+
MODULE PanelProcedures
! ---------------------------------------------------------------------------
! PURPOSE - Define the derived type Panel and collect the procedures that
!  operate on data of type Panel.
USE NetworkProcedures
IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER,PUBLIC:: PANEL_VERSION = '0.51 (27 Dec 2002)'
  
  TYPE:: Panel
    INTEGER:: network  ! index of network to which the panel belongs
    INTEGER:: row,col  ! identify the location within the network
    INTEGER:: compMethod, expMethod
    REAL,DIMENSION(3):: center
    REAL,DIMENSION(3):: normal
    REAL:: area
    REAL:: cosdel   ! radians
    REAL:: cp   ! pressure coefficient
    REAL:: deltar
  END TYPE Panel  
  
  INTEGER,PRIVATE:: DBG=3
  
  PUBLIC::  AverageNetworkProperties  
  PRIVATE:: ComputeCenterAreaNormal
  PUBLIC::  ComputeNetworkForcesAndMoments
  PUBLIC::  ComputePressureCoefficients
  
  PUBLIC::  CountPanels
  
  PUBLIC::  CreatePanels
  PRIVATE:: CrossProduct
  PUBLIC::  DumpPanels
  PUBLIC::  DumpPressures

!----------------------------------------------------------------------------

CONTAINS
!+
SUBROUTINE AverageNetworkProperties(a,b)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the average characteristics
  TYPE(Network),INTENT(IN OUT),DIMENSION(:):: a
  TYPE(Panel),INTENT(IN),DIMENSION(:):: b
  
  INTEGER:: k,n
!----------------------------------------------------------------------------
  DO k=1,SIZE(a)
    a(k)%averageCenter(:)=0.0
    a(k)%averageNormal(:)=0.0
    a(k)%totalArea=0.0
  END DO
  
  DO k=1,SIZE(b)
    n=b(k)%network
    a(n)%averageCenter(:)=a(n)%averageCenter + b(k)%area*b(k)%center(:)
    a(n)%averageNormal(:)=a(n)%averageNormal + b(k)%area*b(k)%normal(:)
    a(n)%totalArea=a(n)%totalArea + b(k)%area
  END DO
  
  DO k=1,SIZE(a)
    IF (a(k)%totalArea <= 0.0) CYCLE
    a(k)%averageCenter(:)=a(k)%averageCenter(:)/a(k)%totalArea
    a(k)%averageNormal(:)=a(k)%averageNormal(:)/a(k)%totalArea
  END DO
    
  RETURN
END Subroutine AverageNetworkProperties   ! ---------------------------------  
 
!+
SUBROUTINE ComputeCenterAreaNormal(x,y,z, i,j, reversed, center,area,normal)
! ---------------------------------------------------------------------------
! PURPOSE - 
  REAL,INTENT(IN),DIMENSION(:,:):: x,y,z
  INTEGER,INTENT(IN):: i,j   ! grid coordinates
  LOGICAL,INTENT(IN):: reversed
  REAL,INTENT(OUT),DIMENSION(3):: center
  REAL,INTENT(OUT):: area
  REAL,INTENT(OUT),DIMENSION(3):: normal

  REAL,DIMENSION(3):: p1,p2,p3,p4,p31,p42

!----------------------------------------------------------------------------
  p1(1)=x(i-1,j-1)
  p1(2)=y(i-1,j-1)
  p1(3)=z(i-1,j-1)
  p2(1)=x(i-1,j)
  p2(2)=y(i-1,j)
  p2(3)=z(i-1,j)
  p3(1)=x(i,j)
  p3(2)=y(i,j)
  p3(3)=z(i,j)
  p4(1)=x(i,j-1)
  p4(2)=y(i,j-1)
  p4(3)=z(i,j-1)
  p31=p3-p1
  p42=p4-p2
  center=0.25*(p1+p2+p3+p4)
  IF (reversed) THEN
    normal=CrossProduct(p42,p31)
  ELSE
    normal=CrossProduct(p31,p42)
  END IF    
  
  area=SQRT(SUM(normal**2)) 
  IF (area==0.0) THEN
    normal=0.0
  ELSE
    normal=normal/area
  END IF
  
  area=0.5*area
  
!  WRITE(DBG,*) 'Panel',i,j
!  WRITE(DBG,'(3F12.5)') p1,p2,p3,p4
!  WRITE(DBG,'(4F15.5)') normal,area
  RETURN
END Subroutine ComputeCenterAreaNormal   ! -----------------------------------  

!+
SUBROUTINE ComputeNetworkForcesAndMoments(nets,panels,xref,yref,zref)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the forces and moments on each network based on the 
!  current pressure coefficients in the panels array.
  TYPE(Network),INTENT(IN OUT),DIMENSION(:):: nets
  TYPE(Panel),INTENT(IN),DIMENSION(:):: panels
  REAL,INTENT(IN):: xref,yref,zref
  
  TYPE(Panel):: a
  INTEGER:: j,k
  REAL:: fmag   ! magnitude of the force vector
  REAL:: xx
!----------------------------------------------------------------------------
  DO k=1,SIZE(nets)
    nets(k)%fx=0.0
    nets(k)%fy=0.0
    nets(k)%fz=0.0
    nets(k)%mx=0.0
    nets(k)%my=0.0
    nets(k)%mz=0.0
  END DO
    
  DO k=1,SIZE(panels)
    a=panels(k)
    j=a%network
    fmag=a%area*a%cp

    xx=-fmag*a%normal(1)
    nets(j)%fx=nets(j)%fx+xx
    nets(j)%my=nets(j)%my+xx*(a%center(3)-zref)
    nets(j)%mz=nets(j)%mz-xx*(a%center(2)-yref)

    xx=-fmag*a%normal(2)
    nets(j)%fy=nets(j)%fy+xx
    nets(j)%mx=nets(j)%mx-xx*(a%center(3)-zref)
    nets(j)%mz=nets(j)%mz+xx*(a%center(1)-xref)

    xx=-fmag*a%normal(3)
    nets(j)%fz=nets(j)%fz+xx
    nets(j)%mx=nets(j)%mx+xx*(a%center(2)-yref)
    nets(j)%my=nets(j)%my-xx*(a%center(1)-xref)
  END DO

  RETURN
END Subroutine ComputeNetworkForcesAndMoments   ! ---------------------------  
  
!+
SUBROUTINE ComputePressureCoefficients(cpstag,cpExpansion,mach,panels)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the pressure coefficient on each panel for a given flight
!  condition.

! PURPOSE - Compute pressure coefficient using various methods of   
!    hypersonic aerodynamics.                                       

USE Methods
  REAL,INTENT(IN):: cpstag   ! stagnation pressure coeff (compMethod=1)
  REAL,INTENT(IN):: cpExpansion
  REAL,INTENT(IN):: mach
  TYPE(Panel),INTENT(IN OUT),DIMENSION(:):: panels
  
  REAL:: cosdel = 0.0
  REAL:: cp
  REAL:: delta ! impact angle, radians
  INTEGER:: k
  REAL,PARAMETER:: PI=3.14159265, HALFPI=0.5*PI

  REAL:: etac=1.0, enpm=1.0, xcent=1.0   ! look in level2.cmn  getpat ??
! ---------------------------------------------------------------------------
  DO k=1,SIZE(panels)
    cosdel=panels(k)%cosdel
    delta=HALFPI-ACOS(cosdel)
    cp=0.0
    IF (delta < 0.0) THEN
      SELECT CASE(panels(k)%expMethod)   ! ----- EXPANSION
        CASE(1)
          cp=0.0
        CASE(2)   ! Newtonian Prandtl Meyer
          cp=NewtonianPrandtlMeyer(mach,delta,cpstag)
        CASE(3)   ! Prandtl-Meyer from free stream   
          cp=PrandtlMeyer(mach,delta)
        CASE(4)
          cp=ConeAtAngleOfAttack()
        CASE(5)
          cp=VanDykeUnified(mach,delta)
        CASE(6) 
          cp=-1.0/mach**2        ! base pressure
        CASE(7)   ! Shock expansion
          cp=0.0  ! later    
        CASE(8)   ! input value
          cp=cpExpansion
        CASE(9)   ! free molecular flow
          cp=0.0  ! later      
        CASE(10)  
          cp=DahlemBuck(mach,delta,cosdel)
        CASE(11)
          cp=ACMempirical(mach,delta)     
        CASE(12)   ! half Prandtl-Meyer from freestream
          cp=0.5*PrandtlMeyer(mach,delta)      
        CASE(13)   ! Tangent Cone (Edwards formulation)
          cp=TangentConeEdwards(mach,delta)
      END SELECT
    ELSE
      SELECT CASE(panels(k)%compMethod)   ! ------ COMPRESSION
        CASE(1)   ! modified Newtonian
          cp=Newtonian(cpstag,delta)
        CASE(2)
          cp=NewtonianPrandtlMeyer(mach,delta,cpstag)
        CASE(3)  
          cp=TangentWedge(mach, delta)
        CASE(4)  
          cp=TangentWedgeInfiniteMach(mach,delta)
        CASE(5)   ! Tangent-Cone Empirical
          cp=OldTangentCone(mach,delta)
        CASE(6)
          cp=ConeAtAngleOfAttack()
        CASE(7)  
          cp=VanDykeUnified(mach,delta)
        CASE(8)
          cp=BluntBodyViscous()
        CASE(9)   ! shock-expansion
          cp=0.0  ! later
        CASE(10)   ! free molecular flow
          cp=0.0   ! later
        CASE(11)   ! use input value of cpstag
          cp=cpstag
        CASE(12)
          cp=HankeyFlatSurface(mach,delta,cosdel)
        CASE(13)
          cp=SmythDeltaWing(mach,delta)
        CASE(14)
          cp=DahlemBuck(mach,delta,cosdel)
        CASE(15)
          cp=BlastWave(mach,cpstag,etac,enpm,xcent)                 
        CASE(16) 
          cp=OSUBluntBody(cpstag,mach,delta)
        CASE(17)
          cp=TangentConeEdwards(mach,delta)
      END SELECT    
    END IF
    panels(k)%cp=cp

  END DO
  RETURN
END Subroutine ComputePressureCoefficients   ! ------------------------------   

!+
FUNCTION CountPanels(a) RESULT(k)   ! a PURE function
! ---------------------------------------------------------------------------
! PURPOSE - Count the number of panels that would be defined be the array of
!  networks called a. Remember that when talking of networks, rows and cols
!  refer to the number of gridpoints. The number of panels in a network of
!  m columns and n rows is (m-1)*(n-1).

  TYPE(Network),INTENT(IN),DIMENSION(:):: a 
 
  INTEGER:: j,k
! ---------------------------------------------------------------------------
  k=0
  DO j=1,SIZE(a)
    k=k+(a(j)%cols-1)*(a(j)%rows-1)
  END DO
  RETURN
END Function CountPanels   ! ------------------------------------------------    
  
!+
SUBROUTINE CreatePanels(networks,panels)
! ---------------------------------------------------------------------------
! PURPOSE - Step thru the array of networks making a panel from each
!  quadrilateral
  TYPE(Network),INTENT(IN),DIMENSION(:):: networks
  TYPE(Panel),INTENT(OUT),DIMENSION(:):: panels
  
  INTEGER:: i,j,k,n
! ---------------------------------------------------------------------------
  k=0
  DO n=1,SIZE(networks)
    DO j=2,networks(n)%cols
      DO i=2,networks(n)%rows
        k=k+1
        panels(k)%network=n
        panels(k)%row=i-1
        panels(k)%col=j-1
        panels(k)%compMethod=networks(n)%compMethod
        panels(k)%expMethod=networks(n)%expMethod
        CALL ComputeCenterAreaNormal(networks(n)%x,networks(n)%y,networks(n)%z, &
          i,j, networks(n)%reversed, panels(k)%center, panels(k)%area, panels(k)%normal)
      END DO
    END DO
  END DO
        
  RETURN
END Subroutine CreatePanels   ! ---------------------------------------------  

!+
FUNCTION CrossProduct(a,b) RESULT(c)
! ---------------------------------------------------------------------------
  REAL,INTENT(IN),DIMENSION(3):: a,b
  
  REAL,DIMENSION(3):: c
!----------------------------------------------------------------------------
  c(1)=a(2)*b(3)-a(3)*b(2)
  c(2)=a(3)*b(1)-a(1)*b(3)
  c(3)=a(1)*b(2)-a(2)*b(1)
  RETURN
END Function CrossProduct   ! -----------------------------------------------

!+
SUBROUTINE DumpPanels(efu, a)
! ---------------------------------------------------------------------------
! PURPOSE - Print the defining data of an array of type Panel
  INTEGER,INTENT(IN):: efu
  TYPE(Panel),INTENT(IN),DIMENSION(:):: a
  
  TYPE(Panel):: b
  CHARACTER(LEN=*),PARAMETER:: FMT = '(2I4, 2I3, 3F10.4, 3F7.3, F12.4)'
  INTEGER:: k
!----------------------------------------------------------------------------
  WRITE(efu,*) 'PANEL DATA'
  WRITE(efu,'(T4,A)') '# net row col   center x,y,x   normal x,y,z   area'
  DO k=1,SIZE(a)
    b=a(k)
    WRITE(efu,FMT) k,b%network,b%row,b%col,b%center, b%normal, b%area
  END DO
  RETURN
END Subroutine DumpPanels   ! -----------------------------------------------

SUBROUTINE DumpPressures(efu,panels)
! ---------------------------------------------------------------------------
  INTEGER,INTENT(IN):: efu
  TYPE(Panel),INTENT(IN),DIMENSION(:):: panels
  
  CHARACTER(LEN=*),PARAMETER:: FMT = '(2I4, 2I3, 3F11.4, 2F9.5, F7.2)'
  INTEGER:: k
  TYPE(Panel):: a  
  REAL,PARAMETER:: PI = 3.14159265
!----------------------------------------------------------------------------
  WRITE(efu,'(T4,A)') &
   '# net  C  E       x          y          z       Cp      cosdel  delta'
  DO k=1,SIZE(panels)
    a=panels(k)
    WRITE(efu,FMT) k,a%network, a%compMethod, a%expMethod, &
      a%center, a%cp, a%cosdel, a%deltar*(180/PI)
  END DO  
  RETURN
END Subroutine DumpPressures   ! --------------------------------------------

END Module PanelProcedures   ! ==============================================

!+
MODULE HypersonicAerodynamicProcedures

IMPLICIT NONE
  INTEGER,PARAMETER,PRIVATE:: DBG = 3

CONTAINS
!+
SUBROUTINE BodyToWind(alpha,beta,phi, cx,cy,cz, cDrag,cyPrime,cLift)
! ---------------------------------------------------------------------------
! PURPOSE - Resolve the aerodynamic forces in body axes in wind axes
IMPLICIT NONE
  REAL,INTENT(IN):: alpha,beta,phi   ! degrees
  REAL,INTENT(IN):: cx,cy,cz
  REAL,INTENT(OUT):: cDrag,cyPrime,cLift
  REAL:: ca,sa, cb,sb, cr,sr
  REAL,PARAMETER:: PI = 3.14159265
!----------------------------------------------------------------------------
  ca=COS((PI/180)*alpha)
  sa=SIN((PI/180)*alpha)
  cb=COS((PI/180)*beta)
  sb=SIN((PI/180)*beta)
  cr=COS((PI/180)*phi)
  sr=SIN((PI/180)*phi)
  
  cDrag = cx*ca*cb - cy*sr*sa*cb -cy*cr*sb + cz*cr*sa*cb - cz*sr*sb
  cLift = -cx*sa-cy*sr*ca+cz*cr*ca
  cyPrime = cx*ca*sb - cy*sr*sa*sb + cy*cr*cb + cz*cr*sa*sb + cz*sr*cb

  RETURN
END Subroutine BodyToWind   ! -----------------------------------------------  

!+
SUBROUTINE ComputeSolutions(networks,panels, cpstag,cpexpansion, &
  mach, sref,cbar,span, xref,yref,zref, solutions)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the pressures on each panel. From these, compute the 
!  forces and moments on each network. Sum these to get the total forces and 
!  moments on the complete configuration at each flight condition
USE NetworkProcedures
USE PanelProcedures
  TYPE(Network),INTENT(IN OUT),DIMENSION(:):: networks
  TYPE(Panel),INTENT(IN OUT),DIMENSION(:):: panels
  REAL,INTENT(IN):: cpstag,cpexpansion
  REAL,INTENT(IN):: mach
  REAL,INTENT(IN):: sref,cbar,span
  REAL,INTENT(IN):: xref,yref,zref
  REAL,INTENT(IN OUT),DIMENSION(:,:):: solutions
!  REAL:: cx,cy,cz, cmx,cmy,cmz

  REAL:: a,b
  CHARACTER(LEN=*),PARAMETER:: FMT = '(A,I3,A,I3,A,F6.2,A,F6.2)'
  INTEGER:: j,k
  INTEGER:: nsol
  REAL,PARAMETER:: PI=3.14159265, HALFPI=0.5*PI
  REAL,DIMENSION(3):: uinf   ! freestream velocity in body coordinate system
!----------------------------------------------------------------------------
  nsol=SIZE(solutions,1)
  DO k=1,nsol
    a=solutions(k,1)*PI/180   ! alpha, radians
    b=solutions(k,2)*PI/180   ! beta, radians
    WRITE(*,FMT)' Computing solution ', k, ' of ', nsol, ' at alpha=', &
         solutions(k,1), '  and beta=', solutions(k,2)
    uinf(1) =  COS(a)*COS(b)
    uinf(2) = -SIN(b)
    uinf(3) =  SIN(a)*COS(b)
    DO j=1,SIZE(panels)
      panels(j)%cosdel = -DOT_PRODUCT(panels(j)%normal, uinf)
      panels(j)%deltar = HALFPI-ACOS(panels(j)%cosdel)
    END DO        

    CALL ComputePressureCoefficients(cpstag,cpExpansion,mach,panels)
    WRITE(DBG,'(A,F9.4,A,F9.4)' ) 'PRESSURE COEFF AT ALPHA=', &
      solutions(k,1), ',   BETA=', solutions(k,2)
    CALL DumpPressures(DBG,panels)
    
    CALL ComputeNetworkForcesAndMoments(networks,panels,xref,yref,zref)
    CALL DumpNetworkForcesAndMoments(DBG,networks)
    solutions(k,3:8)=0.0  
    
    DO j=1,SIZE(networks)
      solutions(k,3)=solutions(k,3)+networks(j)%fx   ! axial force
      solutions(k,4)=solutions(k,4)+networks(j)%fy   ! side force
      solutions(k,5)=solutions(k,5)+networks(j)%fz   ! normal force
      solutions(k,6)=solutions(k,6)+networks(j)%mx   ! rolling moment
      solutions(k,7)=solutions(k,7)+networks(j)%my   ! pitching moment
      solutions(k,8)=solutions(k,8)+networks(j)%mz   ! yawing moment 
    END DO

    CALL BodyToWind(solutions(k,1),solutions(k,2),0.0, &
      solutions(k,3), solutions(k,4), solutions(k,5),  &
      solutions(k,9), solutions(k,10), solutions(k,11) ) ! drag, side', lift
!    cx=cx/sref
!    cy=cy/sref
!    cz=cz/sref
!    cmx=cmx/(span*sref)
!    cmy=cmy/(cbar*sref)
!    cmz=cmz/(span*sref)
!    WRITE(OUT, '(I4,6F12.5)' ) k, cx,cy,cz, cmx,cmy,cmz
!      PrintForcesAndMoments(i, alphaDeg[i], betaDeg[i]);
  END DO
  solutions(:,3:11)=solutions(:,3:11)/sref   ! force coefficients
  solutions(:,6)=solutions(:,6)/span   ! rolling moment coefficient
  solutions(:,7)=solutions(:,7)/cbar   ! pitching moment coefficient
  solutions(:,8)=solutions(:,8)/span   ! yawing moment coefficient
  
  RETURN
END Subroutine ComputeSolutions   !   ---------------------------------------

!+
FUNCTION DetermineNumberOfCases(alpha,beta) RESULT(n)
! ---------------------------------------------------------------------------
! PURPOSE - Determine the number of cases defined by the alpha and beta 
!  arrays. Count backwards from the end of the arrays until an entry is found
!  with either alpha /=0 or beta /= 0. 
  REAL,INTENT(IN),DIMENSION(:):: alpha,beta
  INTEGER:: k,n
!----------------------------------------------------------------------------
  DO k=SIZE(alpha),1,-1
    IF (alpha(k)/=0.0 .OR. beta(k)/=0.0) EXIT
  END DO
  IF (k <= 0) k=1   ! you always want 1, even if all alpha, beta =0
  n=k
  RETURN
END Function DetermineNumberOfCases   ! -------------------------------------  

!+
FUNCTION GetDateTimeStr() RESULT(s)
! ---------------------------------------------------------------------------
! PURPOSE - Return a string with the current date and time
  CHARACTER(LEN=*),PARAMETER:: MONTH="JanFebMarAprMayJunJulAugSepOctNovDec"
  CHARACTER(LEN=*),PARAMETER:: FMT = "(I2.2,A1,I2.2,I3,A3,I4)"
  CHARACTER(LEN=15):: s
  INTEGER,DIMENSION(8):: v
!----------------------------------------------------------------------------
  CALL DATE_AND_TIME(VALUES=v)

  WRITE(s,FMT) v(5), ':', v(6), v(3), MONTH(3*v(2)-2:3*v(2)), v(1)
  RETURN
END FUNCTION GetDateTimeStr   ! ---------------------------------------------

!+
SUBROUTINE PrintMethods(efu)
! ---------------------------------------------------------------------------
! PURPOSE - Print the method codes because the user cannot possibly remember
!   which code goes with each method
  INTEGER,INTENT(IN):: efu   ! external file unit for output
  CHARACTER(LEN=*),PARAMETER:: FMT = '(T4,A)'
!----------------------------------------------------------------------------
  WRITE(efu,*)   'COMPRESSION METHODS'
  WRITE(efu,FMT) ' 1  Modified Newtonian'
  WRITE(efu,FMT) ' 2  Newtonian-Prandtl-Meyer'
  WRITE(efu,FMT) ' 3  Tangent Wedge'
  WRITE(efu,FMT) ' 4  Tangent Wedge Infinite Mach'
  WRITE(efu,FMT) ' 5  Old Tangent Cone'
  WRITE(efu,FMT) ' 6  Cone At Angle Of Attack (later)'
  WRITE(efu,FMT) ' 7  VanDyke Unified'
  WRITE(efu,FMT) ' 8  Blunt Body Viscous (later)'
  WRITE(efu,FMT) ' 9  Shock Expansion (later)'
  WRITE(efu,FMT) '10  Free Molecular Flow (later)'
  WRITE(efu,FMT) '11  Input value of CpStag'
  WRITE(efu,FMT) '12  Hankey Flat Surface'
  WRITE(efu,FMT) '13  Smyth Delta Wing'
  WRITE(efu,FMT) '14  Modified Dahlem-Buck'
  WRITE(efu,FMT) '15  BlastWave (later)'
  WRITE(efu,FMT) '16  OSUBluntBody'
  WRITE(efu,FMT) '17  Tangent Cone (Edwards)'
  WRITE(efu,*)
  WRITE(efu,*) 'EXPANSION METHODS'
  WRITE(efu,FMT) ' 1  Cp=0'
  WRITE(efu,FMT) ' 2  NewtonianPrandtlMeyer'
  WRITE(efu,FMT) ' 3  PrandtlMeyer'
  WRITE(efu,FMT) ' 4  ConeAtAngleOfAttack(later)'
  WRITE(efu,FMT) ' 5  VanDykeUnified'
  WRITE(efu,FMT) ' 6  Vacuum'
  WRITE(efu,FMT) ' 7  Shock Expansion (later)'
  WRITE(efu,FMT) ' 8  Input Value'
  WRITE(efu,FMT) ' 9  Free Molecular Flow(later)'
  WRITE(efu,FMT) '10  Modified Dahlem-Buck'
  WRITE(efu,FMT) '11  ACMempirical(later)'
  WRITE(efu,FMT) '12  half Prandtl-Meyer from freestream'

  RETURN
END Subroutine PrintMethods   ! ---------------------------------------------  


!+
SUBROUTINE PrintSolutions(efu,solutions)
! ---------------------------------------------------------------------------
! PURPOSE - Print the accumulated solutions in body axes and wind axes.
  INTEGER,INTENT(IN):: efu   ! external file unit for output
  REAL,INTENT(IN),DIMENSION(:,:):: solutions
  
  INTEGER:: k,n
!----------------------------------------------------------------------------
  n=SIZE(solutions,1)
  WRITE(efu,*)
  WRITE(efu,*) 'SOLUTIONS FOR COMPLETE VEHICLE IN BODY AXES'  
  WRITE(efu, '(T3,A)' ) &
   '#  alpha   beta      cx        cy        cz       cmx       cmy       cmz'
  WRITE(efu, '(I3,2F7.2,6F10.5)' ) (k, solutions(k,1:8), k=1,n)
  
  WRITE(efu,*)
  WRITE(efu,*) 'SOLUTIONS FOR COMPLETE VEHICLE IN WIND AXES'
  WRITE(efu, '(T3,A)' ) &
   '#  alpha   beta     clift    cdrag     cpitch'
  WRITE(efu, '(I3,2F7.2,3F10.5)' ) (k, solutions(k,1:2), &
    solutions(k,11),solutions(k,9), solutions(k,7), k=1,n)
  WRITE(efu,*)
    
  RETURN
END Subroutine PrintSolutions   ! -------------------------------------------

END Module HypersonicAerodynamicProcedures



PROGRAM HypersonicAerodynamics                             ! \hyper\hyper.f90
! ---------------------------------------------------------------------------
! PURPOSE - Compute aerodynamic forces and moments on an arbitrary vehicle in
!  high speed flow. The local pressure coefficients are computed using various 
!  methods of hypersonic aerodynamics.                                       
! AUTHOR  - Ralph L. Carmichael, Public Domain Aeronautical Software
!  based on original work of Arvel E. Gentry, Douglas N. Smyth and
!  Wayne R. Oliver of Douglas Aircraft.
!  Original work was sponsored by United States Air Force Wright Laboratories
!  Charlotte Craidon, NASA Langley, coordinator of the Langley Wireframe
!                                     Geometry Standard (LaWgs)
         
! REVISION HISTORY                                                  
!   DATE  VERS PERSON  STATEMENT OF CHANGES                         
!   1964?   -    AEG   Original program and inspiration             
!   1973    IV   AEG,DNS,WRO  Release of Mark IV version      
!   1977         CC    Published NASA TM 85767                      
! 30Mar88  0.1   RLC   Original coding (in Pascal)                  
! 24May88  0.2   RLC   Added tangent-wedge empirical                
! 17Mar92  0.3   RLC   Made into a unit for Turbo Pascal            
! 24Mar92  0.31  RLC   Combined constants in one place              
! 12May96  0.4   RLC   Converted to a module for Fortran (90)
! 30Jun96  0.5   RLC   All routines are functions, not subroutines
! 23Nov02  0.6   RLC   Total rewrite in Fortran 95
! 13Dec02  0.65  RLC   Revised input scheme; added more methods
! 31Dec02  0.7   RLC   Last revisions for release in PDAS v.8

USE NetworkProcedures
USE PanelProcedures
USE HypersonicAerodynamicProcedures
IMPLICIT NONE

  CHARACTER(LEN=15):: dateTimeStr
  INTEGER,PARAMETER:: IN=1, WGS=2, DBG=3, OUT=4
  INTEGER,PARAMETER:: MAXCASES = 100, MAXNETS=200
  CHARACTER(LEN=*),PARAMETER:: VERSION = '0.7 (31 Dec 2002)'
  
  REAL,DIMENSION(MAXCASES):: alpha,beta
  REAL:: cbar
  INTEGER,DIMENSION(MAXNETS):: cMethods,eMethods
  REAL:: cpExpansion
  REAL:: cpStag
  CHARACTER(LEN=*),PARAMETER:: FMT1 = '(" Mach=", F6.3,' // &
    ' "  Sref=", F12.4, "  Cbar=", F8.4, "  Span=", F8.4)'
  CHARACTER(LEN=*),PARAMETER:: FMT2 = '(" Moment Reference,   x=",' // &
    ' F12.4, "  y=", F12.4, "  z=", F8.4)'
  CHARACTER(LEN=*),PARAMETER:: FMT3 = '(I6, " networks with ", I7, " panels")'
  CHARACTER(LEN=80):: generalTitle
  TYPE(Network),DIMENSION(MAXNETS):: networks
  TYPE(Panel),ALLOCATABLE,DIMENSION(:):: panels
  INTEGER:: knet
  REAL:: mach
  INTEGER:: nnets,npanels,ncases
  INTEGER,DIMENSION(MAXNETS):: reversed
  REAL,ALLOCATABLE,DIMENSION(:,:):: solutions  ! ncases,11 keep for print
  REAL:: span
  REAL:: sref
  CHARACTER(LEN=80):: title
  CHARACTER(LEN=80):: wgsFileName
  REAL:: xref,yref,zref
  
  NAMELIST /hyp/ wgsFileName, title, cMethods,eMethods,reversed, &
    alpha,beta,cbar,mach,span,sref,xref,yref,zref, cpStag,cpExpansion
!----------------------------------------------------------------------------
  CALL Welcome()
  READ(IN,hyp)
  OPEN(UNIT=WGS,FILE=wgsFileName,STATUS='OLD',ACTION='READ')
  CALL ScanNetworks(WGS, NETS=nnets)
  
  READ(WGS,'(A)') generalTitle
  WRITE(*,*) 'The general title of the WGS file is'
  WRITE(*,*) Trim(generalTitle)
  DO knet=1,nnets
    CALL ReadOneNetwork(WGS,networks(knet))
    networks(knet)%compMethod=cMethods(knet)
    networks(knet)%expMethod=eMethods(knet)
    networks(knet)%reversed=.FALSE.
  END DO  
  DO knet=1,SIZE(reversed)
    IF (reversed(knet)<=0) CYCLE
    IF (reversed(knet)>nnets) CYCLE
    networks(reversed(knet))%reversed=.TRUE.
  END DO  
  WRITE(*,*) nnets, ' networks have been read.'
  CLOSE(UNIT=WGS)
  CALL DumpNetworkGridPoints(DBG, networks(1:nnets))
  
  CALL CreateLocalImageNetworks(nnets,networks)
  WRITE(*,*) 'After local imaging, there are', nnets, ' networks.'
  
  CALL TransformNetworks(networks(1:nnets))
  
  CALL CreateGlobalImageNetworks(nnets,networks)
  WRITE(*,*) 'After global imaging, there are', nnets, ' networks.'
 
  npanels=CountPanels(networks(1:nnets))
  ALLOCATE(panels(npanels)) 
  CALL CreatePanels(networks(1:nnets),panels)
  WRITE(*,*) npanels, ' panels created.'
  CALL DumpPanels(DBG,panels)
  
  CALL AverageNetworkProperties(networks(1:nnets),panels)
  CALL DumpNetworkProperties(DBG, networks(1:nnets))
  
  ncases=DetermineNumberOfCases(alpha,beta)
  ALLOCATE(solutions(ncases,11))
  solutions(1:ncases,1)=alpha(1:ncases)
  solutions(1:ncases,2)=beta(1:ncases)
  
  WRITE(OUT,*) Trim(title)
  WRITE(OUT,FMT1) mach,sref,cbar,span
  WRITE(OUT,FMT2) xref,yref,zref
  WRITE(OUT,FMT3) nnets,npanels
  CALL ComputeSolutions(networks(1:nnets),panels, cpstag,cpexpansion, &
    mach, sref,cbar,span, xref,yref,zref, solutions)

  CALL PrintSolutions(OUT,solutions)
  CALL PrintNetworkMethodsAndNormals(OUT, networks(1:nnets))
  CALL PrintMethods(OUT)
  

  CALL DeAllocateNetworks(networks(1:nnets))
  DEALLOCATE(panels)
  WRITE(*,*) 'Normal termination.'
  STOP
  
CONTAINS


!+
SUBROUTINE Welcome()
! ---------------------------------------------------------------------------
! PURPOSE - Greet user, get name of input file and open it. Open the output
!  file and the log file. Set default values for all input quantities.
USE Methods, ONLY: METHODS_VERSION
  CHARACTER(LEN=*),PARAMETER:: GREETING = "hyper1"
  CHARACTER(LEN=*),PARAMETER:: AUTHOR = &
    "Ralph L. Carmichael, Public Domain Aeronautical Software"
  CHARACTER(LEN=*),PARAMETER:: MODIFIER = "none"

  INTEGER:: errCode
  CHARACTER(LEN=132):: fileName
!----------------------------------------------------------------------------
  WRITE(*,*) GREETING
  WRITE(*,*) "Version "//VERSION                   ! VERSION is global
  WRITE(*,*) AUTHOR
  WRITE(*,*) "Modified by "//Trim(MODIFIER)
  dateTimeStr=GetDateTimeStr()

  DO 
    WRITE(*,*) 'Enter the name of the input file: '
    READ(*,'(A)') fileName
    IF (Len_Trim(fileName)==0) STOP
    OPEN(UNIT=IN, FILE=fileName, STATUS='OLD',   &
      IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
    IF (errCode==0) EXIT
    OPEN(UNIT=IN, FILE=Trim(fileName)//'.inp', STATUS='OLD',   &
      IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
    IF (errCode==0) EXIT
    OPEN(UNIT=IN, FILE=Trim(fileName)//'.dat', STATUS='OLD',   &
      IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
    IF (errCode==0) EXIT
    WRITE(*,*) 'Unable to open this file. Try again.'
  END DO    
  INQUIRE(UNIT=IN, NAME=fileName) 
  WRITE(*,*) "Reading from "//Trim(fileName)
  
  OPEN(UNIT=DBG, FILE='hyper.dbg', STATUS='REPLACE',   &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  IF (errCode == 0) THEN
    WRITE(DBG,*) "--- Created by hyper, version "//VERSION
    WRITE(DBG,*) "--- Uses module NetworkProcedures, version "//NETWORK_VERSION
    WRITE(DBG,*) "--- Uses module PanelProcedures, version "//PANEL_VERSION
    WRITE(DBG,*) "--- Uses module Methods, version "//METHODS_VERSION
    WRITE(DBG,*) "--- Date: "//dateTimeStr
    WRITE(DBG,*) "--- Reading from file: "//Trim(fileName)
  ELSE
    WRITE(*,*) "Unable to create log file. Fatal error"
    STOP
  END IF

  OPEN(UNIT=OUT, FILE='hyper.out', STATUS='REPLACE',   &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  IF (errCode == 0) THEN
    WRITE(OUT,*) "--- Created by hyper, version "//VERSION
    WRITE(OUT,*) "--- Date: "//dateTimeStr
    WRITE(OUT,*) "--- Uses module Methods, version "//METHODS_VERSION
    WRITE(OUT,*) "--- Reading from file: "//Trim(fileName)
  ELSE
    WRITE(*,*) "Unable to create output file. Fatal error"
    STOP
  END IF
  
  xref=0.0   ! set default values for all input data
  yref=0.0
  zref=0.0
  cbar=1.0
  span=1.0
  sref=1.0
  mach=5.0
  alpha(:)=0.0
  beta(:)=0.0
  reversed(:)=0
  cMethods(:)=1
  eMethods(:)=1
  cpStag=2.0
  cpExpansion=0.0

  RETURN
END Subroutine Welcome   ! --------------------------------------------------
  
  
END Program HypersonicAerodynamics   ! ======================================