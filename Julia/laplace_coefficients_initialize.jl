include("laplace_wisdom.jl")

function laplace_coefficients_initialize(jmax::Integer,alpha::Number)
# This computes the Laplace coefficients via recursion.
const nmax=7
b0=Array(eltype(alpha),nmax,jmax+1) # Array to hold the coefficients
# Compute the highest two Laplace coefficients using Wisdom's series approach:
for j=0:jmax
  for i=0:nmax-1
    b0[i+1,j+1]=laplace_wisdom(0.5,j,i,alpha)/alpha^i
  end
end
return b0
end
