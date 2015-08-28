include("laplace_wisdom.jl")

function laplace_coefficients_initialize(jmax::Int64,alpha::Float64)
# This computes the Laplace coefficients via recursion.
nmax=7
b0=zeros(nmax,jmax+1) # Array to hold the coefficients
# Compute the highest two Laplace coefficients using Wisdom's series approach:
for j=0:jmax
  for i=0:nmax-1
    b0[i+1,j+1]=laplace_wisdom(.5,j,i,alpha)/alpha^i
  end
end
b0
#return
end
