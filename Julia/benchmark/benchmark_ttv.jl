data=vec(readdlm("../kepler62ef_planets.txt",',',Float64) ) 
include("../test_ttv.jl")  
@time ttv1,ttv2=test_ttv(5,40,20,data); # This call includes time to compile
Profile.clear_malloc_data()
@time ttv1,ttv2=test_ttv(5,40,20,data,WriteOutput=false,num_evals=100000);

println("# Ignore this: ",hash(ttv1)+hash(ttv2)) # This just makes sure optimizer doesn't optimize away important calculations.

# Test that results haven't changed much
data_inner_out=readdlm("inner_ttv.txt",Float64)
data_outer_out=readdlm("outer_ttv.txt",Float64)
data_inner_ref=readdlm("inner_ttv.txt.ref",Float64)
data_outer_ref=readdlm("outer_ttv.txt.ref",Float64)
max_abs_diff = max(maximum(abs(data_inner_out.-data_inner_ref)),maximum(abs(data_outer_out.-data_outer_ref))) 
println("# Maximum difference: ",max_abs_diff)

