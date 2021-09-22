include("ttv_nplanet.jl")
function ttv_wrapper(tt0,nplanet::Int64,ntrans::Vector{Int64},params::Vector{T},jmax::Integer,EM::Bool) where T<:Real
  # These lines need modification for different choices of parameters:
  if nplanet == 2
    n1,n2 = ntrans
  elseif nplanet == 3
    n1,n2,n3 = ntrans
  elseif nplanet == 4
    n1,n2,n3,n4 = ntrans
  end

  # Call ttv_nplanet:
  ttv = ttv_nplanet(nplanet,jmax,ntrans,params[1:5*nplanet])
  # We measure transit times,not TTVs,so add back in the linear ephemeris:
  t01 = params[3]
  per1 = params[2]
  ttv1 = collect(range(t01,stop = t01+per1*(n1-1),length = n1)) #this doesnt account for skipped transits
  for i=1:n1
    ttv1[i]+= ttv[1,i]
  end
  t02 = params[8]
  per2 = params[7]
  ttv2 = collect(range(t02,stop = t02+per2*(n2-1),length = n2)) #check memory allocation >>>>>>>>>>>>
  for i=1:n2
    if EM
      ttv2[i] += ttv[2,i]
    else
      ts = params[end-2] #tmax sinphi0
      tc = params[end-1] #tmax cosphi0
      deltaphi = params[end]
      ttv2[i] += ttv[2,i] + ts*cos((i-1)*deltaphi) + tc*sin((i-1)*deltaphi)
    end
  end
  # If transit times of additional planets were observable these would need to be added in.
  #println("param2: ",param)
  return [ttv1;ttv2]  #check memory allocation >>>>>>>>>>>>
end

function chisquare(tt0,nplanet,ntrans,params,tt,sigtt,jmax,EM)#,fixp3::Bool = false,p3_cur::Float64 = 0.0)
  chisq = 0.0  #check memory allocation >>>>>>>>>>>>
  # println(params,tt[1],sigtt[1])
  tt_model = ttv_wrapper(tt0,nplanet,ntrans,params,jmax,EM) #,fixp3,p3_cur)
  for j=1:length(tt)
    chisq += (tt[j]-tt_model[j])^2/sigtt[j]^2
  end
  # println(nplanet)
  return chisq
end