# This code calls the first-order eccentric TTV code
# compute_ttv.jl.  Please cite Agol & Deck (2015) if
# you make use of this in your research.

type Planet_plane
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

u(gamma::Float64,c1::Float64,c2::Float64)= ((3+gamma^2)*c1+2*gamma*c2)/(gamma*gamma*(1-gamma*gamma))
# m=+/-1
v(z::Float64,d1::Float64,d2::Float64,m::Int64)= ((m*(1-z*z)+6*z)*d1+(2+z*z)*d2)/(z*(1-z*z)*(z+m)*(z+2*m))

include("compute_ttv.jl")
include("laplace_coefficients_initialize.jl")

function test_ttv(jmax,n1,n2,data)
# Performs a test of the transit_times.jl routine
# Set up planets planar-planet types for the inner and outer planets:
p1=Planet_plane(data[1],data[2],data[3],sqrt(data[4]^2+data[ 5]^2),atan2(data[ 5],data[4]))
p2=Planet_plane(data[6],data[7],data[8],sqrt(data[9]^2+data[10]^2),atan2(data[10],data[9]))
time1 = p1.trans0 + linspace(0,n1-1,n1) * p1.period
time2 = p2.trans0 + linspace(0,n2-1,n2) * p2.period
alpha0=(p1.period/p2.period)^(2//3)
# Initialize the computation of the Laplace coefficients:
b0=laplace_coefficients_initialize(jmax+1,alpha0)
# Define arrays to hold the TTVs:
ttv1=zeros(n1)
ttv2=zeros(n2)
# Define arrays to hold the TTV coefficients and Laplace coefficients:
f1=zeros(jmax+2,5)
f2=zeros(jmax+2,5)
b=zeros(Float64,jmax+2,3)
# Call the compute_ttv code which implements equation (33)
compute_ttv!(jmax,p1,p2,time1,time2,ttv1,ttv2,f1,f2,b,alpha0,b0)
# Write the mean ephemeris and TTV results to two files:
writedlm("inner_ttv.txt",[time1 ttv1])
writedlm("outer_ttv.txt",[time2 ttv2])
ttv1,ttv2
end
