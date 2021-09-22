# Computes TTVs with TTVFaster for N planets with pairwise TTV calculation.
include("compute_ttv.jl")

function ttv_nplanet(nplanet::Int64,jmax::Int64,ntrans::Vector{Int64},data::Vector{T}) where T<:Real
  # Need at least two planets!
  @assert(nplanet>=2)
  # The ntrans vectors should have length nplanet:
  @assert(length(ntrans)==nplanet)
  # Define type of ttv array:
  ttv_el_type = eltype(data) == Float64 ? Float64 : Number
  # Need to create an array to store TTVs with maximum length equal to maximum number
  # of transit times of any planet:
  ntransmax = maximum(ntrans)
  #ttv = zeros(ttv_el_type,nplanet,ntransmax)
  ttv = zeros(T,nplanet,ntransmax) #check memory allocation <<<<<<<<<<<<<<<<<<<
  # Each planet requires 5 elements in data: mass_ratio,period,trans0,ecosw,esinw:
  @assert(length(data)==5*nplanet)
  @assert(jmax>=1)  # Should there be a larger minimum?
  for iplanet=1:nplanet
  # Each planet should have at least 2 transits:
    @assert(ntrans[iplanet]>=2)
  end
  for iplanet=1:nplanet-1
  # The periods of the planets should be ordered from least to greatest:
    if (data[(iplanet-1)*5+2] >= data[iplanet*5+2])
      return ttv
    end
  end
  # Set up planets planar-planet types for all of the planets:
  #planet = Array{Planet_plane_hk}(nplanet)
  #planet = Array{Any}(nplanet)
  # Loop over pairs of planets to compute pairwise TTVs
  # Loop over inner planets:
  #println("Looping over planets in ttv_nplanet:")
  for iplanet=1:nplanet-1
    # Create a Planet_plane_hk type for the inner planet:
    p1=Planet_plane_hk(data[(iplanet-1)*5+1],data[(iplanet-1)*5+2],data[(iplanet-1)*5+3],data[(iplanet-1)*5+4],data[(iplanet-1)*5+5])
    # Create an array of times for the inner planet:
    n1 = ntrans[iplanet]
    time1 = collect(p1.trans0 .+ range(0,stop=n1-1,length=n1) .* p1.period)
    # Loop over outer planets:
    for jplanet=iplanet+1:nplanet
      # Create a Planet_plane_hk type for the outer planet:
      p2=Planet_plane_hk(data[(jplanet-1)*5+1],data[(jplanet-1)*5+2],data[(jplanet-1)*5+3],data[(jplanet-1)*5+4],data[(jplanet-1)*5+5])
      # Create an array of times for the outer planet:
      n2 = ntrans[jplanet]
      time2 = collect(p2.trans0 .+ range(0,stop=n2-1,length=n2) .* p2.period)
      # Define arrays to hold the TTVs:
      ttv1=zeros(T,n1)
      ttv2=zeros(T,n2)
      # Call the compute_ttv code which implements equation (33) from Agol & Deck (2016):
      #    println("Calling compute_ttv")
      compute_ttv!(jmax,p1,p2,time1,time2,ttv1,ttv2)
      #    println("Finished compute_ttv")
      for i=1:n1
        ttv[iplanet,i] += ttv1[i]
      end
      for i=1:n2
        ttv[jplanet,i] += ttv2[i]
      end
    end
  end
  return ttv
end
