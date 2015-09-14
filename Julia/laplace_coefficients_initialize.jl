module LaplaceCoefficients

VERSION < v"0.4-dev" && using Docile

export laplace_coefficients_initialize
export laplace_wisdom
include("laplace_wisdom.jl")

"""
# This computes the Laplace coefficients via recursion.
"""
function initialize(jmax::Integer,alpha::Number)
const nmax=7
b0=Array(eltype(alpha),nmax,jmax+1) # Array to hold the coefficients
# Compute the highest two Laplace coefficients using Wisdom's series approach:
for j=0:jmax
  for i=0:nmax-1
    b0[i+1,j+1]=laplace_wisdom(1//2,j,i,alpha)/alpha^i
  end
end
return b0
end

end # module
