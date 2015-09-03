# Computes TTV coefficients for first-order eccentricity
# solution from Agol & Deck (2015).  Please cite this paper
# if you make use of this in your research.

function ttv_succinct!(jmax::Int64,alpha::Float64,f1::Array{Float64,2},f2::Array{Float64,2},b::Array{Float64,2},alpha0::Float64,b0::Array{Float64,2})

# See simple_solution.pdf 7/16/2015

# Fourth-order Taylor expansion approximation of Laplace coefficients:
dalpha = alpha-alpha0
for i=0:2
  for j=0:jmax
    b[j+1,i+1]=b0[i+1,j+1]+dalpha*(b0[i+2,j+1]+0.5*dalpha*(b0[i+3,j+1]+dalpha/3.0*(b0[i+4,j+1]+dalpha*0.25*b0[i+5,j+1])))
  end
end

sqrtalpha = sqrt(alpha)

# Loop over j:
for j=0:jmax
  # \delta_{j1} (this is indirect coefficient which is only needed for j=1)
  dj1 = j==1 ? 1.0 : 0.0

  # Compute dimensionless frequencies (equation 30):
  beta = j*(1-alpha*sqrtalpha)
  kappa =  beta /  (alpha*sqrtalpha)

  # Compute disturbing function coefficients (equation 31):
  A_j00 = b[j+1,1]
  A_j10 =  alpha* b[j+1,2]
  A_j01 = -(A_j10 + A_j00)
  A_j20 =  alpha*alpha * b[j+1,3]
  A_j11 = -(2*A_j10 + A_j20)
  A_j02 = 2*A_j00 + 4*A_j10 + A_j20
  jd=convert(Float64,j)
  # Inner planet coefficients, in order k=0,-1,1,-2,2 (see Table 1):
  if j >=2
    f1[j+1,1]=alpha*u(beta          ,jd*(    A_j00-alpha*dj1),A_j10-alpha*dj1)
    f1[j+1,2]=alpha*u(beta-1.0      ,jd*(-jd*A_j00-0.5*A_j10+1.5*alpha*dj1),-jd*A_j10-0.5*A_j20+alpha*dj1)
    f1[j+1,3]=alpha*u(beta+1.0      ,jd*( jd*A_j00-0.5*A_j10-0.5*alpha*dj1), jd*A_j10-0.5*A_j20-alpha*dj1)
    f1[j+1,4]=alpha*u(beta-alpha*sqrtalpha,jd*( jd*A_j00-0.5*A_j01-2.0*alpha*dj1), jd*A_j10-0.5*A_j11-2.0*alpha*dj1)
    f1[j+1,5]=alpha*u(beta+alpha*sqrtalpha,jd*(-jd*A_j00-0.5*A_j01),-jd*A_j10-0.5*A_j11)
  else
    if j==0
      f1[j+1,4]=alpha*u(beta-alpha*sqrt(alpha),jd*( jd*A_j00-0.5*A_j01-2.0*alpha*dj1), jd*A_j10-0.5*A_j11-2.0*alpha*dj1)
    else
      f1[j+1,1]=alpha*u(beta          ,jd*(    A_j00-alpha*dj1),A_j10-alpha*dj1)
      f1[j+1,2]=alpha*u(beta-1.0      ,jd*(-jd*A_j00-0.5*A_j10+1.5*alpha*dj1),-jd*A_j10-0.5*A_j20+alpha*dj1)
      f1[j+1,3]=alpha*u(beta+1.0      ,jd*( jd*A_j00-0.5*A_j10-0.5*alpha*dj1), jd*A_j10-0.5*A_j20-alpha*dj1)
      f1[j+1,4]=alpha*u(beta-alpha*sqrt(alpha),jd*( jd*A_j00-0.5*A_j01-2.0*alpha*dj1), jd*A_j10-0.5*A_j11-2.0*alpha*dj1)
    end
  end
  # Add in the k=\pm 1 coefficients (note that d1 & d2 are the same as c1 & c2 for k=0):
  if j >= 1
    f1[j+1,2]=f1[j+1,2]+alpha*v(beta,jd*(A_j00-alpha*dj1),A_j10-alpha*dj1,-1)
    f1[j+1,3]=f1[j+1,3]+alpha*v(beta,jd*(A_j00-alpha*dj1),A_j10-alpha*dj1, 1)
  end

# Now for the outer planet:
# Outer planet coefficients, in order k=0,-2,2,-1,1 (see Table 1):
# TODO: Test if it is worth making an one_over_alpha_squared?
  if j >= 2
    f2[j+1,1]=u(kappa,-jd*(A_j00-dj1/alpha^2),A_j01-dj1/alpha^2)
    f2[j+1,2]=u(kappa-1,-jd*(jd*A_j00-0.5*A_j01-0.5*dj1/alpha^2),jd*A_j01-0.5*A_j02-dj1/alpha^2)
    f2[j+1,3]=u(kappa+1,-jd*(-jd*A_j00-0.5*A_j01+1.5*dj1/alpha^2),-jd*A_j01-0.5*A_j02+dj1/alpha^2)
    f2[j+1,4]=u(kappa-1/(alpha*sqrtalpha),-jd*(-jd*A_j00-0.5*A_j10),-jd*A_j01-0.5*A_j11)
  else
    if j == 1
      f2[j+1,1]=u(kappa,-jd*(A_j00-dj1/alpha^2),A_j01-dj1/alpha^2)
      f2[j+1,2]=u(kappa-1,-jd*(jd*A_j00-0.5*A_j01-0.5*dj1/alpha^2),jd*A_j01-0.5*A_j02-dj1/alpha^2)
      f2[j+1,3]=u(kappa+1,-jd*(-jd*A_j00-0.5*A_j01+1.5*dj1/alpha^2),-jd*A_j01-0.5*A_j02+dj1/alpha^2)
    end
  end
  f2[j+1,5]=u(kappa+1/(alpha*sqrtalpha),-jd*(jd*A_j00-0.5*A_j10-2.0*dj1/alpha^2),jd*A_j01-0.5*A_j11-2.0*dj1/alpha^2)
# Add in the k=\pm 2 coefficients (note that d1 & d2 are the same as c1 & c2 for k=0):
  if j >= 1
    f2[j+1,2]=f2[j+1,2]+v(kappa,-jd*(A_j00-dj1/alpha^2),A_j01-dj1/alpha^2,-1)
    f2[j+1,3]=f2[j+1,3]+v(kappa,-jd*(A_j00-dj1/alpha^2),A_j01-dj1/alpha^2, 1)
  end
# That's it!
end
return 
end
