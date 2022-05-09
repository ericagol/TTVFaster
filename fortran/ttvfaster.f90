MODULE ttvfaster
	
	! Fortran90 version of TTVFaster (Kat Deck & Eric Agol) by David Kipping
	! Please cite the paper (Agol & Deck 2015) if you make use of this code in
	! your research.

 implicit none

 CONTAINS

	! Try to preserve original for loop indices, but update memory allocations +1

 ! =======================================================
 SUBROUTINE compute_ttv(nt1,nt2,time1,time2,param,jmax,model1,model2)

 implicit none
 
 ! Inputs 
 INTEGER, INTENT(IN) :: nt1, nt2, jmax
 REAL(8), DIMENSION(nt1), INTENT(IN) :: time1
 REAL(8), DIMENSION(nt2), INTENT(IN) :: time2
 REAL(8), DIMENSION(10), INTENT(IN) :: param
 ! Intermediate
 INTEGER :: i, j
 REAL(8) :: p1massratio, p1trans0, p1period, p1ecos, p1esin
 REAL(8) :: p2massratio, p2trans0, p2period, p2ecos, p2esin
 REAL(8) :: alpha
 REAL(8), DIMENSION(nt1) :: lam11, lam21, psi1
 REAL(8), DIMENSION(nt2) :: lam12, lam22, psi2
 REAL(8), DIMENSION(5,jmax+2) :: f1, f2
 REAL(8), DIMENSION(nt1) :: ttv1
 REAL(8), DIMENSION(nt2) :: ttv2
 ! Output
 REAL(8), DIMENSION(nt1), INTENT(OUT) :: model1
 REAL(8), DIMENSION(nt2), INTENT(OUT) :: model2
 ! Constants
 REAL(8), PARAMETER :: twopi = 6.283185307179586D0
 
 ! Planet 1 parameters
 p1massratio = param(1)
 p1trans0 = param(2)
 p1period = param(3)
 p1ecos = param(4)
 p1esin = param(5)
 ! Planet 2 parameters
 p2massratio = param(6)
 p2trans0 = param(7)
 p2period = param(8)
 p2ecos = param(9)
 p2esin = param(10)
 
 ! Compute the semi-major axis ratio of the planets
 alpha = DABS(p1period/p2period)**(2.0/3.0)
 ! Compute the longitudes of the planets at times of transit of planet 1 (equation 32)
 DO i=1,nt1
	 lam11(i) = twopi*(time1(i)-p1trans0)/p1period + 2.0D0*p1esin
	 lam21(i) = twopi*(time1(i)-p2trans0)/p2period + 2.0D0*p2esin
	 ! Compute the difference in in longitudes at times of transit of planet 1
	 psi1(i) = lam11(i) - lam21(i)
 END DO
 ! Compute the longitudes of the planets at times of transit of planet 2 (equation 32)
 DO i=1,nt2
	 lam12(i) = twopi*(time2(i)-p1trans0)/p1period + 2.0D0*p1esin
	 lam22(i) = twopi*(time2(i)-p2trans0)/p2period + 2.0D0*p2esin
	 ! Compute the difference in in longitudes at times of transit of planet 2
	 psi2(i) = lam12(i) - lam22(i)
 END DO
 ! Compute the coefficients (need one higher than jmax)
 call ttv_succint(jmax+1,alpha,f1,f2)
 
 ttv1(:) = 0.0D0
 ttv2(:) = 0.0D0
 ! Compute TTVs for inner planet (equation 33)
 DO j=1,jmax ! python does range(1,jmax+1) = 1,2,...,jmax
	 ttv1(:) = ttv1(:) + f1(1+0,1+j)*DSIN(j*psi1(:))
	 ttv1(:) = ttv1(:) + f1(1+1,1+j)*( p1ecos*DSIN( (j-1)*lam11(:) - j*lam21(:) ) + &
	                                   p1esin*DCOS( (j-1)*lam11(:) - j*lam21(:) ) )
	 ttv1(:) = ttv1(:) + f1(1+2,1+j)*( p1ecos*DSIN( (j+1)*lam11(:) - j*lam21(:) ) - &
																     p1esin*DCOS( (j+1)*lam11(:) - j*lam21(:) ) )									 
	 ttv1(:) = ttv1(:) + f1(1+3,1+j-1)*( p2ecos*DSIN( (j-1)*lam11(:) - j*lam21(:) ) + &
															  	     p2esin*DCOS( (j-1)*lam11(:) - j*lam21(:) ) )
	 ttv1(:) = ttv1(:) + f1(1+4,1+j+1)*( p2ecos*DSIN( (j+1)*lam11(:) - j*lam21(:) ) - &
																       p2esin*DCOS( (j+1)*lam11(:) - j*lam21(:) ) )
 END DO
 ! Now multiply by the mass ratio and divide by mean motion
 ttv1(:) = ttv1(:)*p1period*p2massratio/twopi

 ! Compute TTVs for outer planet (equation 33)
 DO j=1,jmax
	 ttv2(:) = ttv2(:) + f2(1+0,1+j)*DSIN(j*psi2(:))
	 ttv2(:) = ttv2(:) + f2(1+1,1+j)*( p2ecos*DSIN( j*lam12(:) - (j+1)*lam22(:) ) + &
	                                   p2esin*DCOS( j*lam12(:) - (j+1)*lam22(:) ) )
	 ttv2(:) = ttv2(:) + f2(1+2,1+j)*( p2ecos*DSIN( j*lam12(:) - (j-1)*lam22(:) ) - &
																     p2esin*DCOS( j*lam12(:) - (j-1)*lam22(:) ) )									 
	 ttv2(:) = ttv2(:) + f2(1+3,1+j+1)*( p1ecos*DSIN( j*lam12(:) - (j+1)*lam22(:) ) + &
															  	     p1esin*DCOS( j*lam12(:) - (j+1)*lam22(:) ) )
	 ttv2(:) = ttv2(:) + f2(1+4,1+j-1)*( p1ecos*DSIN( j*lam12(:) - (j-1)*lam22(:) ) - &
																       p1esin*DCOS( j*lam12(:) - (j-1)*lam22(:) ) )
 END DO
 ! Now multiply by the mass ratio and divide by mean motion
 ttv2(:) = ttv2(:)*p2period*p1massratio/twopi
 
 ! Add the TTVs to the ephemeris and return to user
 DO i=1,nt1
   model1(i) = p1trans0 + (i-1)*p1period + ttv1(i)
 END DO
 DO i=1,nt2
   model2(i) = p2trans0 + (i-1)*p2period + ttv2(i)
 END DO
 												 
 END SUBROUTINE compute_ttv
 ! =======================================================

 ! =======================================================
 FUNCTION ufn(gammz,z1,z2)
 ! Define function u (equation 34)
	 
 implicit none

 REAL(8), INTENT(IN) :: gammz, z1, z2
 REAL(8) :: ufn
 
 ufn = ( ( (3.0+gammz**2)*z1 + 2.0*gammz*z2 ) / &
       gammz**2 / (1.0-gammz**2) )
  
 RETURN

 END FUNCTION
 ! =======================================================
 
 ! =======================================================
 FUNCTION vfn(zeta,d1,d2,m)
 ! Define function v (equation 34)
	 
 implicit none

 INTEGER, INTENT(IN) :: m ! m = +/- 1
 REAL(8), INTENT(IN) :: zeta, d1, d2
 REAL(8) :: vfn
 
 vfn = ((m*(1.0 - zeta**2) + 6.0*zeta)*d1 + & 
       (2.0 + zeta**2)*d2) / (zeta*(1.0-zeta**2)* &
	     (zeta+m)*(zeta+2.0*m))
  
 RETURN

 END FUNCTION
 ! =======================================================

 ! =======================================================
 FUNCTION laplace_wisdom(s,i,j,a)

 ! Code due to Jack Wisdom.
 ! Compute Laplace coefficients and Leverrier derivatives
 !       j
 !  j   d     i
 ! a   ---   b (a)
 !       j    s
 !     da
 ! by series summation.

 REAL(8), INTENT(IN) :: s, a
 INTEGER, INTENT(IN) :: j
 INTEGER :: i
 REAL(8), PARAMETER :: laplaceeps = 1.0D-12
 REAL(8) :: asquare, term
 REAL(8) :: factorone, factortwo, factorthree, factorfour
 INTEGER :: q0, q, k
 REAL(8) :: laplace_wisdom

 asquare = a*a
 i = ABS(i)
 
 ! compute first term in sum
 IF( j .LE. i ) THEN
   factorfour = 1.0
   DO k=0,j-1 ! python was range(j), so 0,1,2,..,j-1
     factorfour = factorfour*DBLE(i-k)
   END DO
   laplace_wisdom = factorfour
   q0 = 0
 ELSE
   q0 = (j+1-i)/2
   laplace_wisdom = 0.0
   factorfour = 1.0	
 END IF
 
 ! compute factors for terms in sum
 factorone = s
 factortwo = s + DBLE(i)
 factorthree = DBLE(i) + 1.0
 q = 1
 ! no contribution from q = 0
 DO q=1,q0-1 ! python is range (1,q0) so 1,2..,q0-1
   factorone   = factorone*(   s + DBLE(q)     )
   factortwo   = factortwo*(   s + DBLE(i+q)   )
   factorthree = factorthree*(     DBLE(i+q+1) )
 END DO
 
 term = asquare*factorone*factortwo/(factorthree*DBLE(q))
 
 ! sum series
 DO WHILE (term*factorfour .GT. laplaceeps)
   factorfour = 1.0
   DO k=0,j-1 ! python was range (j), so 0,1,2,..,j-1
     factorfour = factorfour*( 2.0*DBLE(q) + DBLE(i-k) )
   END DO
   laplace_wisdom = laplace_wisdom + term*factorfour
   factorone   = factorone   + 1.0
   factortwo   = factortwo   + 1.0
   factorthree = factorthree + 1.0
   q = q + 1
   term = term*asquare*factorone*factortwo/(factorthree*DBLE(q))
 END DO
 
 ! fix coefficient
 DO k=0,i-1 ! python was range(i) so 0,1,2,..,i-1
   laplace_wisdom = laplace_wisdom*(s+DBLE(k)) / (DBLE(k)+1.0)
 END DO
 
 IF( q0 .LE. 0 ) THEN
   laplace_wisdom = laplace_wisdom*2.0*a**i
 ELSE
   laplace_wisdom = laplace_wisdom*2.0*a**(2*q0 + i - 2)
 END IF
 
 RETURN

 END FUNCTION
 ! =======================================================
 
 ! =======================================================
 SUBROUTINE ttv_succint(jmax,alpha,f1,f2)

 implicit none
 
 INTEGER :: j, k
 INTEGER, INTENT(IN) :: jmax
 REAL(8), INTENT(IN) :: alpha
 REAL(8) :: Pratio, beta, kappa, dj1
 REAL(8) :: Ajaa, Ajba, Ajab, Ajca, Ajbb, Ajac
 REAL(8), DIMENSION(5) :: c1, c2, gammy
 REAL(8), DIMENSION(3,jmax+1) :: b
 REAL(8), DIMENSION(5,jmax+1), INTENT(OUT) :: f1, f2

 ! Succinct form for coefficients of first-order TTV formula
 !    Parameters
 !    ----------
 !    jmax:  maximum value of j over which to compute the coefficients
 !    alpha:  a_1/a_2 of the two planets
 !    Returns
 !    -------
 !    b:  Matrix of Laplace coefficients.
 !    f1:  Coefficients for the inner planet.  For each value of
 !        j=0 to jmax, there are 5 coefficients:
 !        the f_{1,j}^(0), f_{1,j}^(+-1), and
 !        f_{1,j}^(+-2) coefficients.
 !        The +-1 coefficients correspond to
 !        j*(lam_1-lam_2)+-(lam_1-omega_1) arguments.
 !        The +-2 coefficients correspond to
 !        j*(lam_1-lam_2)+-(lam_1-omega_2) arguments.
 !    f2:  Coefficients for the outer planet.  For each value of
 !        j=0 to jmax, there are 5 coefficients:
 !        the f_{1,j}^(0), f_{1,j}^(+-2), and
 !        f_{1,j}^(+-1) coefficients.
 !        The +-1 coefficients correspond to
 !        j*(lam_1-lam_2)+-(lam_2-omega_1) arguments.
 !        The +-2 coefficients correspond to
 !        j*(lam_1-lam_2)+-(lam_2-omega_2) arguments.

 f1(:,:) = 0.0D0
 f2(:,:) = 0.0D0
 Pratio = DSQRT(alpha**3)

 ! Compute the Laplace coefficients
 call laplace_coefficients(jmax,alpha,b)
 
 ! Now loop over j values
 DO j=0,jmax ! python uses range(jmax+1)=0,1,2,...,jmax
   IF( j .EQ. 1 ) THEN
	   dj1 = 1.0
   ELSE
	   dj1 = 0.0
   END IF
   beta = j*( 1.0 - Pratio )
   kappa = beta/Pratio
   
   ! Compute the disturbing function [note we have to increase memory indices by 1, and switch order]
   Ajaa = b(1+0,1+j)
   Ajba = alpha*b(1+1,1+j)
   Ajab = -(Ajba + Ajaa)
   Ajca = alpha**2*b(1+2,1+j)
   Ajbb = -(2.0*Ajba + Ajca)
   Ajac = 2.0*Ajaa + 4.0*Ajba + Ajca
 
   ! Inner planer coefficients
   gammy(1) = beta
   gammy(2) = beta - 1.0
   gammy(3) = beta + 1.0
   gammy(4) = beta - Pratio
   gammy(5) = beta + Pratio
   
   c1(1) = alpha*j*( Ajaa*1.0 - 0.5*Ajab*0.0 - 0.5*Ajba*0.0 - alpha*dj1*1.0 )
   c1(2) = alpha*j*(-Ajaa*j   - 0.5*Ajab*0.0 - 0.5*Ajba*1.0 + alpha*dj1*1.5 )
   c1(3) = alpha*j*( Ajaa*j   - 0.5*Ajab*0.0 - 0.5*Ajba*1.0 - alpha*dj1*0.5 )
   c1(4) = alpha*j*( Ajaa*j   - 0.5*Ajab*1.0 - 0.5*Ajba*0.0 - alpha*dj1*2.0 )
   c1(5) = alpha*j*(-Ajaa*j   - 0.5*Ajab*1.0 - 0.5*Ajba*0.0 - alpha*dj1*0.0 )
   
   c2(1) = alpha*( Ajba*1.0 - 0.5*Ajbb*0.0 - 0.5*Ajca*0.0 - alpha*dj1*1.0 )
   c2(2) = alpha*(-Ajba*j   - 0.5*Ajbb*0.0 - 0.5*Ajca*1.0 + alpha*dj1*1.0 )
   c2(3) = alpha*( Ajba*j   - 0.5*Ajbb*0.0 - 0.5*Ajca*1.0 - alpha*dj1*1.0 )
   c2(4) = alpha*( Ajba*j   - 0.5*Ajbb*1.0 - 0.5*Ajca*0.0 - alpha*dj1*2.0 )
   c2(5) = alpha*(-Ajba*j   - 0.5*Ajbb*1.0 - 0.5*Ajca*0.0 - alpha*dj1*0.0 )
   
   IF( j .GE. 2 ) THEN
     DO k=0,4 ! python does range(5) = 0,1,2,3,4
		   f1(1+k,1+j) = ufn( gammy(1+k), c1(1+k), c2(1+k) )
	   END DO
   ELSE
	   IF( j .EQ. 0 ) THEN
		   f1(1+3,1+j) = ufn( gammy(1+3), c1(1+3), c2(1+3) )
	   ELSE
	     DO k=0,3 ! python does range(4) = 0,1,2,3
		     f1(1+k,1+j) = ufn( gammy(1+k), c1(1+k), c2(1+k) )
	     END DO  
     END IF
	 END IF
	 
   ! Add in the k=\pm 1 coefficients (note that d1 & d2 are the same as)
   ! c1 & c2 for k=0
   IF( j .GE. 1 ) THEN
	   f1(1+1,1+j) = f1(1+1,1+j) + vfn( beta, c1(1+0), c2(1+0),-1 )
	   f1(1+2,1+j) = f1(1+2,1+j) + vfn( beta, c1(1+0), c2(1+0), 1 )
   END IF
   
   ! Now for the outer planet
   ! Outer planet coefficients, in order k=0,-2,2,-1,1 (see Table 1)
	 ! [manual assignments]
	 gammy(1) = kappa
	 gammy(2) = kappa - 1.0
	 gammy(3) = kappa + 1.0
	 gammy(4) = kappa - 1.0/Pratio
	 gammy(5) = kappa + 1.0/Pratio
	 
	 c1(1) = -j*( Ajaa*1.0 - 0.5*Ajab*0.0 - 0.5*Ajba*0.0 - dj1/alpha**2*1.0 )
	 c1(2) = -j*( Ajaa*j   - 0.5*Ajab*1.0 - 0.5*Ajba*0.0 - dj1/alpha**2*0.5 ) 
	 c1(3) = -j*(-Ajaa*j   - 0.5*Ajab*1.0 - 0.5*Ajba*0.0 + dj1/alpha**2*1.5 ) 
	 c1(4) = -j*(-Ajaa*j   - 0.5*Ajab*0.0 - 0.5*Ajba*1.0 - dj1/alpha**2*0.0 ) 
	 c1(5) = -j*( Ajaa*j   - 0.5*Ajab*0.0 - 0.5*Ajba*1.0 - dj1/alpha**2*2.0 )
	 
	 c2(1) = ( Ajab*1.0 - 0.5*Ajbb*0.0 - 0.5*Ajac*0.0 - dj1/alpha**2*1.0 )
	 c2(2) = ( Ajab*j   - 0.5*Ajbb*0.0 - 0.5*Ajac*1.0 - dj1/alpha**2*1.0 ) 
	 c2(3) = (-Ajab*j   - 0.5*Ajbb*0.0 - 0.5*Ajac*1.0 + dj1/alpha**2*1.0 )
	 c2(4) = (-Ajab*j   - 0.5*Ajbb*1.0 - 0.5*Ajac*0.0 - dj1/alpha**2*0.0 )
	 c2(5) = ( Ajab*j   - 0.5*Ajbb*1.0 - 0.5*Ajac*0.0 - dj1/alpha**2*2.0 )
	 
	 IF( j .GE. 2 ) THEN
		 DO k=0,3 ! python does range(4) = 0,1,2,3
			 f2(1+k,1+j) = ufn( gammy(1+k), c1(1+k), c2(1+k) )
		 END DO
	 ELSE
		 IF( j .EQ. 1 ) THEN
			 DO k=0,2 ! python does range(3) = 0,1,2
				 f2(1+k,1+j) = ufn( gammy(1+k), c1(1+k), c2(1+k) )
			 END DO
		 END IF
	 END IF
	 f2(1+4,1+j) = ufn(gammy(1+4),c1(1+4),c2(1+4))
	 ! Add in the k=\pm 2 coefficients (note that d1 & d2 are the same as
   ! c1 & c2 for k=0)
	 IF( j .GE. 1 ) THEN
	   f2(1+1,1+j) = f2(1+1,1+j) + vfn( kappa, c1(1+0), c2(1+0),-1)
		 f2(1+2,1+j) = f2(1+2,1+j) + vfn( kappa, c1(1+0), c2(1+0), 1)
	 END IF
 END DO
 
 END SUBROUTINE ttv_succint
 ! =======================================================
 
 ! =======================================================
 SUBROUTINE laplace_coefficients(j,alpha,b)
	 
 !This computes the Laplace coefficients via recursion.
 ! j spans 0 to jmax
 ! b(1,j+1), b(2,j+1), b(3,j+1) go out

 implicit none
 
 INTEGER, INTENT(IN) :: j
 REAL(8), INTENT(IN) :: alpha
 INTEGER :: j0
 REAL(8) :: jd
 REAL(8), DIMENSION(3,j+1), INTENT(OUT) :: b

 ! Because of fortran indexing, we have to just add +1 onto
 ! the indices uses by the b array. We swap the indices around
 ! to optimize memory usage

 ! Compute the highest two Laplace coefficients using
 ! Wisdom's series approach
 b(:,:) = 0.0D0
 b(1+0,1+j)   = laplace_wisdom(0.5D0,j,0,alpha)
 b(1+0,1+j-1) = laplace_wisdom(0.5D0,j-1,0,alpha)
 b(1+1,1+j)   = laplace_wisdom(0.5D0,j,1,alpha)/alpha
 b(1+1,1+j-1) = laplace_wisdom(0.5D0,j-1,1,alpha)/alpha
 b(1+2,1+j)   = laplace_wisdom(0.5D0,j,2,alpha)/alpha/alpha
 b(1+2,1+j-1) = laplace_wisdom(0.5D0,j-1,2,alpha)/alpha/alpha
 
 ! The rest can be found with the recurion formula
 j0 = j - 2
 DO WHILE (j0 .GE. 0)
   ! Recurrence relations (derived on 9/9/14; also see Brouwer & Clemence)
   jd = DBLE(j0)
   b(1+0,1+j0) = (-(2.0*jd+3.0)*b(1+0,1+j0+2) + &
                 (1.0+alpha*alpha)/alpha*2.0*(jd+1.0)*b(1+0,1+j0+1))/(2.0*jd+1.0)
   ! Recurrence for first derivative
   b(1+1,1+j0) = b(1+1,1+j0+2) - (2.0*alpha*(jd+1.0)*b(1+0,1+j0+1) - &
                 jd*b(1+0,1+j0) - (jd+2.0)*b(1+0,1+j0+2) )/alpha
   ! Recurrence for second derivative	
   b(1+2,1+j0) = b(1+2,1+j0+2) - 2.0*(jd+1.0)*b(1+1,1+j0+1) - &
                 ( jd*b(1+0,1+j0) + &
			           (jd+2.0)*b(1+0,1+j0+2) )/alpha/alpha + &
			           ( jd*b(1+1,1+j0) + (jd+2.0)*b(1+1,1+j0+2) )/alpha
   j0 = j0 - 1
 END DO
 
 END SUBROUTINE laplace_coefficients
 ! =======================================================

END MODULE ttvfaster
