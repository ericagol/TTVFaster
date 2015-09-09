# Computes transit timing variations to linear order
# in eccentricity.  Please cite Agol & Deck (2015) if
# you make use of this in published research.

#include("ttv_coeff.jl")
include("ttv_succinct.jl")

immutable Planet_plane
# Parameters of a planet in a plane-parallel system
# Mass ratio of the planet to the star:
  mass_ratio :: Float64
# Initial time of transit:
  period   :: Float64
  trans0   :: Float64
  eccen    :: Float64
# longitude of periastron measured from line of sight, in radians:
  omega    :: Float64
end

u(gamma::Float64,c1::Float64,c2::Float64)= ((3+gamma*gamma)*c1+2*gamma*c2)/(gamma*gamma*(1-gamma*gamma))
# m=+/-1
v(z::Float64,d1::Float64,d2::Float64,m::Int64)= ((m*(1-z*z)+6*z)*d1+(2+z*z)*d2)/(z*(1-z*z)*(z+m)*(z+2*m))

function compute_ttv!(jmax::Int64,p1::Planet_plane,p2::Planet_plane,time1::Array{Float64,1},time2::Array{Float64,1},ttv1::Array{Float64,1},ttv2::Array{Float64,1},f1::Array{Float64,2},f2::Array{Float64,2},b::Array{Float64,2},alpha0::Float64,b0::Array{Float64,2})

# Computes transit-timing variations to linear order in
# eccentricity for non-resonant, plane-parallel planets.
# Input:
#   jmax:  Maximum j over which to sum the TTV calculation for both planets
#     p1:  Planet type for inner planet
#     p2:  Planet type for outer planet
#  time1:  Transit times for inner planet
#  time2:  Transit times for outer planet
# alpha0:  Initial alpha for Taylor expansion of coefficients
#     b0:  Laplace coefficients and derivatives for use in Taylor expansion
# Output:
#   ttv1: TTVs of the inner planet
#   ttv2: TTVs of the outer planet
#     f1: TTV coefficients for inner planet
#     f2: TTV coefficients for outer planet
#      b: Laplace coefficients (& derivatives) for outer planet
# Compute the semi-major axis ratio of the planets:
# println(p1.period,p2.period)
const alpha = (p1.period/p2.period)^(2//3)  # Julia supports rational numbers!
# Number of times:
const ntime1 = length(time1)
const ntime2 = length(time2)
# Compute the coefficients:
ttv_succinct!(jmax+1,alpha,f1,f2,b,alpha0,b0)  # I need to compute coefficients one higher than jmax
# Compute TTVs for inner planet (equation 33):
# Compute since of \pomegas:
sin1om=sin(p1.omega)
sin2om=sin(p2.omega)
cos1om=cos(p1.omega)
cos2om=cos(p2.omega)
# Compute mean motions:
n1=2pi/p1.period
n2=2pi/p2.period
# Compute initial longitudes:
lam10=-n1*p1.trans0 + 2*p1.eccen*sin1om
lam20=-n2*p2.trans0 + 2*p2.eccen*sin2om
@inbounds for i=1:ntime1
# Compute the longitudes of the planets at times of transit of planet 1 (equation 49):
  lam11 = n1*time1[i]+lam10
  lam21 = n2*time1[i]+lam20
  psi1  = lam11-lam21 # Compute difference in longitudes at times of transit of planet 1
  sinpsi1=sin(psi1)
  cospsi1=cos(psi1)
  sinlam1om1=sin(lam11-p1.omega)
  coslam1om1=cos(lam11-p1.omega)
  sinlam1om2=sin(lam11-p2.omega)
  coslam1om2=cos(lam11-p2.omega)
  ttv1[i]=0.0
  sinjm1psi1=0.0
  cosjm1psi1=1.0
# Sum over j:
  for j=1:jmax
    sinjpsi1=sinjm1psi1*cospsi1+cosjm1psi1*sinpsi1
    cosjpsi1=cosjm1psi1*cospsi1-sinjm1psi1*sinpsi1
    ttv1[i] += f1[j+1,1]*sinjpsi1
    ttv1[i] += f1[j+1,2]*p1.eccen*(sinjpsi1*coslam1om1-cosjpsi1*sinlam1om1)
    ttv1[i] += f1[j+1,3]*p1.eccen*(sinjpsi1*coslam1om1+cosjpsi1*sinlam1om1)
    ttv1[i] += f1[j  ,4]*p2.eccen*(sinjpsi1*coslam1om2-cosjpsi1*sinlam1om2)
    ttv1[i] += f1[j+2,5]*p2.eccen*(sinjpsi1*coslam1om2+cosjpsi1*sinlam1om2)
    sinjm1psi1=sinjpsi1
    cosjm1psi1=cosjpsi1
  end
# Multiply by period and mass ratio, and divide by 2*Pi:
  ttv1[i] = ttv1[i]*p1.period*p2.mass_ratio/(2pi)
end
# Compute TTVs for outer planet (equation 33):
@inbounds for i=1:ntime2
# Compute the longitudes of the planets at times of transit of planet 2:
  lam12 = n1*time2[i]+lam10
  lam22 = n2*time2[i]+lam20
  psi2  = lam12-lam22 # Compute difference in longitudes at times of transit of planet 2
  sinpsi2=sin(psi2)
  cospsi2=cos(psi2)
  sinlam2om1=sin(lam22-p1.omega)
  coslam2om1=cos(lam22-p1.omega)
  sinlam2om2=sin(lam22-p2.omega)
  coslam2om2=cos(lam22-p2.omega)
  ttv2[i]=0.0
  sinjm1psi2=0.0
  cosjm1psi2=1.0
# Sum over j:
  for j=1:jmax
    sinjpsi2=sinjm1psi2*cospsi2+cosjm1psi2*sinpsi2
    cosjpsi2=cosjm1psi2*cospsi2-sinjm1psi2*sinpsi2
    ttv2[i] += f2[j+1,1]*sinjpsi2
    ttv2[i] += f2[j+1,2]*p2.eccen*(sinjpsi2*coslam2om2-cosjpsi2*sinlam2om2)
    ttv2[i] += f2[j+1,3]*p2.eccen*(sinjpsi2*coslam2om2+cosjpsi2*sinlam2om2)
    ttv2[i] += f2[j+2,4]*p1.eccen*(sinjpsi2*coslam2om1-cosjpsi2*sinlam2om1)
    ttv2[i] += f2[j  ,5]*p1.eccen*(sinjpsi2*coslam2om1+cosjpsi2*sinlam2om1)
    sinjm1psi2=sinjpsi2
    cosjm1psi2=cosjpsi2
  end
# Multiply by period and mass ratio, and divide by 2*Pi:
  ttv2[i] = ttv2[i]*p2.period*p1.mass_ratio/(2pi)
end
# Finished!
return
end

