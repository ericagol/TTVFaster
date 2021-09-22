"""
    Transit_Struct
Structure to hold arrays and other quantities for computing transit. 
# Members
- `r::T`: radius ratio
- `b::T`: impact parameter
- `jmax::Int64`: maximum number of terms in series expansions of ``I_v`` and ``J_v``

"""

mutable struct Transit_Struct{T<:Real}
	tt0::Array{T,1}
	ttv::Array{T,1}
	jmax::Int64
end

n = Transit_Struct()