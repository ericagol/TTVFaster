function laplace_wisdom,s,i,j,a
; double laplace(double s, int i, int j, double a);


;/* Code due to Jack Wisdom */
;/* compute Laplace coefficients and Leverrier derivatives
;          j
;     j   d     i
;    a   ---   b (a)
;          j    s
;        da
;
;   by series summation */

;#define LAPLACE_EPS 1.0e-12
LAPLACE_EPS=1.0d-12

as = a*a

;if (i lt 0) then i = -i
i=abs(i)

if (j le i) then begin     ;/* compute first term in sum */
  factor4 = 1d0
  for k=0, j-1 do factor4 *= double(i - k)
  sum = factor4
  q0 = 0
endif else begin
  q0 = (j + 1 - i)/2
  sum = 0d0
  factor4 = 1d0
endelse

;  /* compute factors for terms in sum */

factor1 = s
factor2 = s + double(i)
factor3 = double(i) + 1d0
for q=1,q0-1 do begin ;/* no contribution for q = 0 */
  factor1 *= s + double(q)
  factor2 *= s + double(i + q)
  factor3 *= double(i+q+1)
endfor
term = as * factor1 * factor2 / (factor3 * double(q))

;  /* sum series */

while ((term*factor4) gt LAPLACE_EPS) do begin
  factor4 = 1d0
  for k=0,j-1 do factor4 *= (2d0*double(q) + double(i - k))
  sum += term * factor4
  factor1 += 1d0
  factor2 += 1d0
  factor3 += 1d0
  q = q+1
  term *= as * factor1 * factor2 / (factor3 * double(q))
endwhile

;  /* fix coefficient */

for k=0,i-1 do sum *= (s + double(k))/(double(k)+1d0)

if(q0 le 0) then sum *= 2d0 * a^i else sum *= 2d0 * a^(2*q0 + i - 2)
return,sum
end
