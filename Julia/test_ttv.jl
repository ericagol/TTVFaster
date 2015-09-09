# This code calls the first-order eccentric TTV code
# compute_ttv.jl.  Please cite Agol & Deck (2015) if
# you make use of this in your research.

include("compute_ttv.jl")
include("laplace_coefficients_initialize.jl")

function test_ttv(jmax::Integer,n1::Integer,n2::Integer,data::Array{Float64,2}; WriteOutput::Bool = true, num_evals::Integer = 1)
@assert(jmax>=1)  # Should there be a larger minimum?
@assert(n1>2)
@assert(n2>2)
@assert(length(data)==10)
# Performs a test of the transit_times.jl routine
# Set up planets planar-planet types for the inner and outer planets:
p1=Planet_plane(data[1],data[2],data[3],sqrt(data[4]^2+data[ 5]^2),atan2(data[ 5],data[4]))
p2=Planet_plane(data[6],data[7],data[8],sqrt(data[9]^2+data[10]^2),atan2(data[10],data[9]))
time1 = collect(p1.trans0 + linspace(0,n1-1,n1) * p1.period)
time2 = collect(p2.trans0 + linspace(0,n2-1,n2) * p2.period)
alpha0=(p1.period/p2.period)^(2//3)
# Initialize the computation of the Laplace coefficients:
b0=laplace_coefficients_initialize(jmax+1,alpha0)
# Define arrays to hold the TTVs:
ttv1=Array(Float64,n1)
ttv2=Array(Float64,n2)
# Define arrays to hold the TTV coefficients and Laplace coefficients:
f1=Array(Float64,jmax+2,5)
f2=Array(Float64,jmax+2,5)
b=Array(Float64,jmax+2,3)
for i in 1:num_evals
   # Call the compute_ttv code which implements equation (33)
   compute_ttv!(jmax,p1,p2,time1,time2,ttv1,ttv2,f1,f2,b,alpha0,b0)
end
if WriteOutput
   # Write the mean ephemeris and TTV results to two files:
   writedlm("inner_ttv.txt",[time1 ttv1])
   writedlm("outer_ttv.txt",[time2 ttv2])
end
   return ttv1,ttv2
end
