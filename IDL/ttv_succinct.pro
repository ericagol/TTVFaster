; This routine computes the coefficients for the first-order
; transit timing variations in Agol & Deck (2015).  Please
; cite the paper if you make use of this code in your research.

; Define u & v functions (equation 34):
function u,gamma,c1,c2
return,((3d0+gamma^2)*c1+2d0*gamma*c2)/gamma^2/(1d0-gamma^2)
end

function v,zeta,d1,d2,m
; m=+/-1
return,((m*(1d0-zeta^2)+6d0*zeta)*d1+(2d0+zeta^2)*d2)/(zeta*(1d0-zeta^2)*(zeta+m)*(zeta+2d0*m))
end

pro ttv_succinct,jmax,alpha,f1,f2,b
; Succinct form for coefficients of first-order TTV formula
; Input:
;   jmax:  maximum value of j over which to compute the coefficients
;  alpha:  a_1/a_2 of the two planets
; Output:
;      b:  Matrix of Laplace coefficients.
;     f1:  Coefficients for the inner planet.  For each value of j=0 to jmax,
;          there are 5 coefficients: the f_{1,j}^(0), f_{1,j}^(+-1), and f_{1,j}^(+-2) coefficients.
;          The +-1 coefficients correspond to j*(lam_1-lam_2)+-(lam_1-omega_1) arguments.
;          The +-2 coefficients correspond to j*(lam_1-lam_2)+-(lam_1-omega_2) arguments.
;     f2:  Coefficients for the outer planet.  For each value of j=0 to jmax,
;          there are 5 coefficients: the f_{1,j}^(0), f_{1,j}^(+-2), and f_{1,j}^(+-1) coefficients
;          The +-1 coefficients correspond to j*(lam_1-lam_2)+-(lam_2-omega_1) arguments.
;          The +-2 coefficients correspond to j*(lam_1-lam_2)+-(lam_2-omega_2) arguments.

f1=dblarr(jmax+1,5)
f2=dblarr(jmax+1,5)

; Compute the Laplace coefficients:
laplace_coefficients3,jmax, alpha,b

; Now loop over j values:
for j=0,jmax do begin
  if(j eq 1) then dj1 = 1d0 else dj1 = 0d0
  beta = j*(1-alpha^1.5)
  kappa =  beta /  alpha^1.5

; Compute the disturbing function coefficients A_jmn (equatio 31):
  A_j00 = b[j,0]
  A_j10 =  alpha* b[j,1]
  A_j01 = -(A_j10 + A_j00)
  A_j20 =  alpha^2 * b[j,2]
  A_j11 = -(2*A_j10 + A_j20)
  A_j02 = 2*A_j00 + 4*A_j10 + A_j20

; Inner planet coefficients, in order k=0,-1,1,-2,2 (see Table 1):
  gamma = beta+[0,-1,1,-alpha^1.5,alpha^1.5]
  c1=alpha*j*(A_j00*[1,-j,j,j,-j]-.5*A_j01*[0,0,0,1,1]-.5*A_j10*[0,1,1,0,0]-alpha*dj1*[1,-1.5,.5,2,0])
  c2=alpha*  (A_j10*[1,-j,j,j,-j]-.5*A_j11*[0,0,0,1,1]-.5*A_j20*[0,1,1,0,0]-alpha*dj1*[1,-1, 1,2,0])
  if(j ge 2) then for k=0,4 do f1[j,k]=u(gamma[k],c1[k],c2[k]) else $
  if(j eq 0) then f1[j,3]=u(gamma[3],c1[3],c2[3]) else begin
    for k=0,3 do f1[j,k]=u(gamma[k],c1[k],c2[k])
  endelse
; Add in the k=\pm 1 coefficients (note that d1 & d2 are the same as c1 & c2 for k=0):
  if(j ge 1) then f1[j,1]=f1[j,1]+v(beta,c1[0],c2[0],-1)
  if(j ge 1) then f1[j,2]=f1[j,2]+v(beta,c1[0],c2[0], 1)
; Now for the outer planet:
; Outer planet coefficients, in order k=0,-2,2,-1,1 (see Table 1):
  gamma = kappa+[0,-1,1,-1d0/alpha^1.5,1d0/alpha^1.5]
  c1=-j*(A_j00*[1,j,-j,-j,j]-.5*A_j01*[0,1,1,0,0]-.5*A_j10*[0,0,0,1,1]-dj1/alpha^2*[1,.5,-1.5,0,2])
  c2=   (A_j01*[1,j,-j,-j,j]-.5*A_j11*[0,0,0,1,1]-.5*A_j02*[0,1,1,0,0]-dj1/alpha^2*[1,1,-1,0,2])
  if(j ge 2) then for k=0,3 do f2[j,k]=u(gamma[k],c1[k],c2[k]) else $
  if(j eq 1) then for k=0,2 do f2[j,k]=u(gamma[k],c1[k],c2[k])
  f2[j,4]=u(gamma[4],c1[4],c2[4])
; Add in the k=\pm 2 coefficients (note that d1 & d2 are the same as c1 & c2 for k=0):
  if(j ge 1) then f2[j,1]=f2[j,1]+v(kappa,c1[0],c2[0],-1)
  if(j ge 1) then f2[j,2]=f2[j,2]+v(kappa,c1[0],c2[0], 1)
; That's it!
endfor
return
end
