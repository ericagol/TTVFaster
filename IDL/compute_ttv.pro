; This routine computes the first-order
; transit timing variations in Agol & Deck (2015).  Please
; cite the paper if you make use of this code in your research.

pro compute_ttv,time,param,model
; Computes transit-timing variations to linear order in
; eccentricity for non-resonant, plane-parallel planets
; to first order in mass ratio and eccentricity.  It is in a format
; that is compatible with optimization with mpcurvefit.pro (Levenberg-Marquardt
; solver).
;
; Input:
;    jmax:  maximum j to evalued - it is passed as a common block variable.
;   param:  contains the (mass ratio,t0,period,e*cos(omega),e*sin(omega)) for each planet
; Output:
;   model:  Contains the transit timing model with the first set for the inner
;           planet and second set for the outer planet.
common order,jmax
; Insert the parameters of the model into a structure for each planet:
p1 = {mass_ratio:param[0],trans0:param[1],period:param[2],eccen:sqrt(param[3]^2+param[4]^2),omega:atan(param[4],param[3])}
p2 = {mass_ratio:param[5],trans0:param[6],period:param[7],eccen:sqrt(param[8]^2+param[9]^2),omega:atan(param[9],param[8])}
; Count the number of times to be computed:
ntime = n_elements(time)
; The inner planet is assumed to be the first set of times
; and the outer planet the second set of times.  The time
; which is less than the prior time, then, is the starting
; time of the outer planet:
indx=where(time[0:ntime-2] gt time[1:ntime-1])
; Make arrays of the transit times for each planet:
time1=time[0:indx] & time2=time[indx+1:ntime-1]
;
twopi = 2d0*!dpi
; Compute the semi-major axis ratio of the planets:
alpha = (p1.period/p2.period)^(2./3.)
; Compute the longitudes of the planets at times of transit of planet 1 (equation 32):
lam11 = twopi*(time1-p1.trans0)/p1.period + 2*p1.eccen*sin(p1.omega)
lam21 = twopi*(time1-p2.trans0)/p2.period + 2*p2.eccen*sin(p2.omega)
;print,lam21
; Compute the longitudes of the planets at times of transit of planet 2 (equation 32):
lam12 = twopi*(time2-p1.trans0)/p1.period + 2*p1.eccen*sin(p1.omega)
lam22 = twopi*(time2-p2.trans0)/p2.period + 2*p2.eccen*sin(p2.omega)
psi1  = lam11-lam21 ; Compute difference in longitudes at times of transit of planet 1
psi2  = lam12-lam22 ; Compute difference in longitudes at times of transit of planet 2
; Number of times:
ntime1 = n_elements(time1)
ntime2 = n_elements(time2)
; Compute the coefficients:
ttv_succinct,jmax+1,alpha,f1,f2,b  ; I need to compute coefficients one higher than jmax
; Create two arrays to hold the TTVs with zeros:
ttv1 = dblarr(ntime1)
ttv2 = dblarr(ntime2)

; Compute TTVs for inner planet (equation 33):
for j=1,jmax do begin
  ttv1 += f1[j,0]*sin(j*psi1)
  ttv1 += f1[  j,1]*p1.eccen*sin((j-1)*lam11-j*lam21+p1.omega)
  ttv1 += f1[  j,2]*p1.eccen*sin((j+1)*lam11-j*lam21-p1.omega)
  ttv1 += f1[j-1,3]*p2.eccen*sin((j-1)*lam11-j*lam21+p2.omega)
  ttv1 += f1[j+1,4]*p2.eccen*sin((j+1)*lam11-j*lam21-p2.omega)
endfor
; Now multiply by the mass ratio and divide by mean motion:
ttv1 = ttv1*p1.period*p2.mass_ratio/(twopi)
; Compute TTVs for outer planet (equation 33):
for j=1,jmax do begin
  ttv2 += f2[j,0]*sin(j*psi2)
  ttv2 += f2[  j,1]*p2.eccen*sin(j*lam12-(j+1)*lam22+p2.omega)
  ttv2 += f2[  j,2]*p2.eccen*sin(j*lam12-(j-1)*lam22-p2.omega)
  ttv2 += f2[j+1,3]*p1.eccen*sin(j*lam12-(j+1)*lam22+p1.omega)
  ttv2 += f2[j-1,4]*p1.eccen*sin(j*lam12-(j-1)*lam22-p1.omega)
endfor
; Now multiply by the mass ratio and divide by mean motion:
ttv2 = ttv2*p2.period*p1.mass_ratio/(twopi)

; Compute secular terms from the appendix (right now this is not used):
e1 = p1.eccen
omega1 = p1.omega
theta1 = lam11 + 2.*e1*sin(lam11-p1.omega)
dtheta1_dt = twopi/p1.period*(1. + 2.*e1*cos(lam11-omega1))
e2 = p2.eccen
omega2 = p2.omega
theta2 = lam22 + 2.*e2*sin(lam22-omega2)
dtheta2_dt = twopi/p2.period*(1. + 2.*e2*cos(lam22-omega2))

A_011   =  - alpha *(2*b[0,1] +  alpha *b[0,2])
A_000   =  b[0,0]
A_100   =  b[1,0]
A_010   =  alpha*b[0,1]
A_110   =  alpha*b[1,1]
A_020   =  alpha^2*b[0,2]
A_120   =  alpha^2*b[1,2]
A_001   =  -(A_010 + A_000)
A_101   =  -(A_110 + A_100)
A_002   =  2.*A_000 + 4.*A_010 + A_020
A_111   = -(2*A_110+A_120)

; Add in equations 64 & 68:
re_c1  = .5*(3.*A_010+A_020)*e1*cos(lam11-omega1) + (2.*A_100+    A_101+A_110+.5*A_111)*e2*cos(lam11-omega2)  ; This is Re[c_1 e^{i\lambda_1}]
im_c1p = .5*(2.*A_010+A_020)*e1*sin(lam11-omega1) + (   A_100+0.5*A_101+A_110+.5*A_111)*e2*sin(lam11-omega2)  ; This is Im[c_1 e^{i\lambda_1}]
dtheta1_sec = p2.mass_ratio*alpha*(-lam11*.25*A_010 - lam11*re_c1 + im_c1p)
dt1_sec = -dtheta1_sec/dtheta1_dt
re_c2  =  .5*(3.*A_001+A_002)*e2*cos(lam22-omega2) + (2.*A_100+    A_110+A_101+.5*A_111)*e1*cos(lam22-omega1)
im_c2p =  .5*(2.*A_001+A_002)*e2*sin(lam22-omega2) + (   A_100+0.5*A_110+A_101+.5*A_111)*e1*sin(lam22-omega1)
dtheta2_sec = p1.mass_ratio*(-lam22*.25*A_001 - lam22*re_c2 + im_c2p)
dt2_sec = -dtheta2_sec/dtheta2_dt

; Add the TTVs to the ephemeris and return to user:
model=[p1.trans0+dindgen(ntime1)*p1.period+ttv1,p2.trans0+dindgen(ntime2)*p2.period+ttv2]
; The following line can be used if the secular terms are desired:
;model=[p1.trans0+dindgen(ntime1)*p1.period+ttv1+dt1_sec,p2.trans0+dindgen(ntime2)*p2.period+ttv2+dt2_sec]
; The following line can be used if only the TTVs are desired:
;model=[ttv1,ttv2]
; Finished!
return
end
