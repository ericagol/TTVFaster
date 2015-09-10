#module TargetDist

param_ref=vec(readdlm("../kepler62ef_planets.txt",',',Float64) ) 
include("../test_ttv.jl")  
calc_ttvs_given_param(param::Array) = vcat(test_ttv(5,40,20,param,WriteOutput=false,num_evals=1));

#=
function set_target_data_obs()
  global const ttv_ref = calc_ttvs_given_param(param_ref)
  global const sigma_ttv = 10/(60*24)*ones(length(ttv_ref))
  param_sgn = sign(param_ref)
  global const param_pert = param_sgn.*exp( log(abs(param_ref)).+0.001*randn(length(param_ref)) )
  nothing
end
set_target_data_obs()
=#

ttv_ref = calc_ttvs_given_param(param_ref)
sigma_ttv = 10/(60*24)*ones(length(ttv_ref))
param_sgn = sign(param_ref)
param_pert = param_sgn.*exp( log(abs(param_ref)).+0.001*randn(length(param_ref)) )


function calc_chisq(theta::Array) 
   ttv_pred = calc_ttvs_given_param(theta)
   chisq = zero(eltype(theta))
   for i in 1:length(ttv_pred)
     chisq += ((ttv_pred[i]-ttv_ref[i]::Float64)/sigma_ttv[i]::Float64)^2
   end
   return chisq
end

log_likelihood(theta::Array) = -0.5*calc_chisq(theta)

using ForwardDiff
my_cache = ForwardDiffCache() # make new cache to pass in to our function

# (grad, allresults) = ForwardDiff.gradient(log_likelihood,param_ref, AllResults, cache=my_cache);
log_target(theta::Array) = log_likelihood(theta)
grad_log_target(theta::Array) = ForwardDiff.gradient(log_likelihood,theta, cache=my_cache);

using Lora
mcmodel = model(log_target,grad=grad_log_target, init=param_pert)

mcchain = run(mcmodel, HMC(0.75), SerialMC(nsteps=1000, burnin=100))




function test_forw_diff(param_ref::Vector)
@time chisq = calc_chisq(param_ref); 
@time chisq_list = typeof(chisq)[calc_chisq(param_pert) for i in 1:100]; 
println("# Ignore this (Chisq):", hash(chisq_list))

@time grad = ForwardDiff.gradient(log_likelihood, param_ref, cache=my_cache); 
@time grad_list = typeof(grad)[ForwardDiff.gradient(log_likelihood, param_pert, cache=my_cache) for i in 1:100]; 
println("# Ignore this (grad): ", hash(grad_list))

@time (grad, allresults) = ForwardDiff.gradient(log_likelihood,param_ref, AllResults, cache=my_cache);
@time grad_list = typeof((grad,allresults))[ForwardDiff.gradient(log_likelihood,param_pert, AllResults, cache=my_cache) for i in 1:100];
println("# Ignore this (grad): ", hash(grad_list))

#@time jac = ForwardDiff.jacobian(log_likelihood,param_ref); 
#@time jac_list = typeof(jac)[ForwardDiff.jacobian(calc_chisq,param_pert) for i in 1:100];
##println("# Ignore this (jac): ", hash(jac_list))

#@time hes_func = ForwardDiff.hessian(log_likelihood);
#@time hes = hes_func(param_ref); 
#@time hes_list = typeof(hes)[hes_func(param_pert) for i in 1:100];
#println("# Ignore this: ", hash(hes_list))

#@time hes, allresults = ForwardDiff.hessian(log_likelihood,param_ref,AllResults); 
#@time hes, allresults = ForwardDiff.hessian(log_likelihood,param_pert,AllResults); 
#println("# Ignore this (hes): ", hash(hesc_list))

return grad_list
end # func test_forw_diff

#end # module TargetDist

# output = test_forw_diff(param_ref)

