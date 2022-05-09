PROGRAM example
use ttvfaster

 implicit none
 
 ! Constants
 INTEGER, PARAMETER :: nt1 = 40
 INTEGER, PARAMETER :: nt2 = 20
 INTEGER, PARAMETER :: sdim = 10
 INTEGER, PARAMETER :: jmax = 6
 ! Inputs 
 REAL(8), DIMENSION(sdim) :: datap, param
 REAL(8), DIMENSION(nt1) :: time1, model1
 REAL(8), DIMENSION(nt2) :: time2, model2
 INTEGER :: i

 ! read-in system parameters
 open(30,FILE='kepler62ef_planets.txt',FORM='FORMATTED',STATUS='UNKNOWN')
 DO i=1,sdim
   read(30,*) datap(i)
 END DO
 close(30)
 
 ! Setup a uniform ephemeris to compute TTVs
 DO i=1,nt1
   time1(i) = datap(2) + (i-1)*datap(3)
 END DO
 DO i=1,nt2
   time2(i) = datap(7) + (i-1)*datap(8)
 END DO
 
 ! Re-arrange planet parameers for ingest to compute_ttv
 param(:) = datap(:)
 
 ! model contains the transit time
 call compute_ttv(nt1,nt2,time1,time2,param,jmax,model1,model2)
 
 ! Save the result
 open(31,FILE='inner_ttv.txt',FORM='FORMATTED',STATUS='UNKNOWN')
 DO i=1,nt1
   write(31,*) model1(i)
 END DO
 close(31)
 open(32,FILE='outer_ttv.txt',FORM='FORMATTED',STATUS='UNKNOWN')
 DO i=1,nt2
   write(32,*) model2(i)
 END DO
 close(32)

END PROGRAM example
