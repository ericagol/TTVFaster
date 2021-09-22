pro laplace_coefficients3,j,alpha,b
; This computes the Laplace coefficients via recursion.
;na=n_elements(alpha)
b=dblarr(j+1,3) ; I'm going to compute coefficient dependence on b^(j) & b^(j-1), as well as derivatives:
; Compute the highest two Laplace coefficients using Wisdom's series approach:
b[j  ,0]=laplace_wisdom(.5d0,j  ,0,alpha)
b[j-1,0]=laplace_wisdom(.5d0,j-1,0,alpha)
b[j,  1]=laplace_wisdom(.5d0,j  ,1,alpha)/alpha
b[j-1,1]=laplace_wisdom(.5d0,j-1,1,alpha)/alpha
b[j,  2]=laplace_wisdom(.5d0,j  ,2,alpha)/alpha^2
b[j-1,2]=laplace_wisdom(.5d0,j-1,2,alpha)/alpha^2

; The rest can be found with the recursion formulae:
j0=j-2
while(j0 ge 0) do begin
; Recurrence relations (derived on 9/9/14; also see Brouwer & Clemence):
  jd = double(j0)
  b[j0,0]=(-(2d0*jd+3d0)*b[j0+2,0]+(1d0+alpha^2)/alpha*2d0*(jd+1d0)*b[j0+1,0])/(2d0*jd+1d0)
; Recurrence for first derivative:
  b[j0,1]=b[j0+2,1]-(2d0*alpha*(jd+1d0)*b[j0+1,0]-jd*b[j0,0]-(jd+2d0)*b[j0+2,0])/alpha
; Recurrence for second derivative:
  b[j0,2]=b[j0+2,2]-2d0*(jd+1d0)*b[j0+1,1]-(jd*b[j0,0]+(jd+2d0)*b[j0+2,0])/alpha^2+$
           (jd*b[j0,1]+(jd+2d0)*b[j0+2,1])/alpha
  j0=j0-1
endwhile
return
end
