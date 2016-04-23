This code was written by Eric Agol 09-07-2015
agol@uw.edu
and edited/optimized by Eric Ford.

Here is an example of using the code in Julia:

```
ulia$ julia  
               _
   _       _ _(_)_     |  A fresh approach to technical computing  
  (_)     | (_) (_)    |  Documentation: http://docs.julialang.org  
   _ _   _| |_  __ _   |  Type "help()" for help.  
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 0.3.0 (2014-08-20 20:43 UTC)  
 _/ |\__'_|_|_|\__'_|  |  Official http://julialang.org/ release  
|__/                   |  x86_64-apple-darwin13.3.0  

julia> data=readdlm("kepler62ef_planets.txt",',',Float64)
1x10 Array{Float64,2}: 3.02306e-5 122.386 -16.5926 -0.00127324 0.0026446 1.67874e-5 267.307 155.466 -0.0025544 0.00117917

julia> include("test_ttv.jl")
test_ttv (generic function with 1 method)

julia> @time ttv1,ttv2=test_ttv(5,40,20,data);
elapsed time: 0.345652398 seconds (13941760 bytes allocated)
julia> @time ttv1,ttv2=test_ttv(5,40,20,data);
elapsed time: 0.000526126 seconds (20404 bytes allocated)
```

This computes the TTVs for a system similar to Kepler-62e/f stored in the file kepler62ef_planet.txt. The TTVs will be written to the files inner_ttv.txt and outer_ttv.txt, as well as stored in the variables ttv1 and ttv2. The test_ttv.jl routine accepts jmax (the maximum j to sum to, in this example 5), ntime1 (number of transits of the inner planet), ntime2 (the number of transits of the outer planet), and data which contains the parameters of the planet.

The file kepler62ef_planets.txt contains a comma-separated set of 10 parameters that describe the system: \mu,Period,t0, e cos(omega), e sin(omega) for each planet, where \mu is the mass ratio of the planet to the star, t0 is the initial transit time (of the averaged orbit), Period is the mean orbital period, e is the eccentricity, and omega is the longitude of periastron.
