param_ref  = vec(readdlm("../kepler62ef_planets.txt",',',Float64) ) 
include("../test_ttv.jl")

function calc_ttvs_given_param(param::Vector) 
  output = test_ttv(5,40,20,param,WriteOutput=false,num_evals=1)
  return vcat( output[1], output[2] )
end

function set_target_data_obs(param::Vector)
  global const ttv_ref = calc_ttvs_given_param(param)
  global const sigma_ttv = 10/(60*24)*ones(length(ttv_ref))
  nothing
end

set_target_data_obs(param_ref)

param_pert_scale = 0.01
param_pert = sign(param_ref).*exp( log(abs(param_ref)).+param_pert_scale*randn(length(param_ref)) )

#=
ttv_ref    = calc_ttvs_given_param(param_ref)
sigma_ttv  = 10/(60*24)*ones(length(ttv_ref))
param_sgn  = sign(param_ref)
param_pert = param_sgn.*exp( log(abs(param_ref)).+0.001*randn(length(param_ref)) )
=#

function calc_chisq(theta::Vector)
   ttv_pred = calc_ttvs_given_param(theta)
   @assert( size(ttv_pred) == size(ttv_ref) )
   chisq = zero(eltype(theta))
   for i in 1:length(ttv_pred)
     delta = (ttv_pred[i]-ttv_ref[i])/sigma_ttv[i]
     chisq += delta*delta
   end
   return chisq
end

log_likelihood(theta::Vector) = -0.5*calc_chisq(theta)

using ForwardDiff
my_cache = ForwardDiff.ForwardDiffCache() # make new cache to pass in to our function

# (grad, allresults) = ForwardDiff.gradient(log_likelihood,param_ref, AllResults, cache=my_cache);
# log_target(theta::Vector) = log_likelihood(theta)
# grad_log_target(theta::Vector) = ForwardDiff.gradient(log_likelihood,theta, cache=my_cache);

function test_forw_diff(param_ref::Vector{Float64}; test_hessian::Bool = false)

println("# Testing evaluating function")
@time ttvs = calc_ttvs_given_param(param_ref); 
@time ttvs_list = typeof(ttvs)[calc_ttvs_given_param(param_pert) for i in 1:1000]; 
println("# Ignore this (ttvs):", hash(ttvs_list))
return

println("# Testing evaluating function")
@time chisq = calc_chisq(param_ref); 
@time chisq_list = typeof(chisq)[calc_chisq(param_pert) for i in 1:1000]; 
println("# Ignore this (Chisq):", hash(chisq_list))

println("# Testing gradient")
@time grad = ForwardDiff.gradient(log_likelihood, param_ref, cache=my_cache); 
@time grad_list = typeof(grad)[ForwardDiff.gradient(log_likelihood, param_pert, cache=my_cache) for i in 1:100]; 
println("# Ignore this (gradient): ", hash(grad_list))

@time (grad, allresults) = ForwardDiff.gradient(log_likelihood,param_ref, AllResults, cache=my_cache);
@time grad_results_list = typeof((grad,allresults))[ForwardDiff.gradient(log_likelihood,param_pert, AllResults, cache=my_cache) for i in 1:100];
println("# Ignore this (eval & gradient): ", hash(grad_results_list))

println("# Testing Jacobbian")
@time jac = ForwardDiff.jacobian(calc_ttvs_given_param,param_ref, cache=my_cache);
@time jac_list = typeof(jac)[ForwardDiff.jacobian(calc_ttvs_given_param,param_ref, cache=my_cache) for i in 1:100];
println("# Ignore this (jacobian): ", hash(jac_list))

if test_hessian
println("# Testing Hessian")
@time hes_func = ForwardDiff.hessian(log_likelihood);
@time hes = hes_func(param_ref); 
@time hes_list = typeof(hes)[hes_func(param_pert) for i in 1:100];
println("# Ignore this (hessian): ", hash(hes_list))

println("# Testing everything at once")
@time hes, allresults = ForwardDiff.hessian(log_likelihood,param_ref,AllResults); 
@time hes, allresults = ForwardDiff.hessian(log_likelihood,param_pert,AllResults); 
println("# Ignore this (hessian, all at once): ", hash(hesc_list))
end
return true
end # func test_forw_diff

test_forw_diff(param_ref)
# test_forw_diff(param_ref, test_hessian=true);  # Doesn't work (yet)

