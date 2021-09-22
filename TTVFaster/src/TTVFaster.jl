"""
    TTVFaster

Computes first order eccentricity transit timing variations (TTVs) with respect to the following initial conditions: the planet-star mass ratio [μ], the initial transit time (of the averaged orbit) [t0], the mean orbital period [Per], the eccentricity [e], and the longitude of periastron [barred ω].
"""

module TTVFaster

# include("ttv_nplanet.jl")
include("ttv_wrapper.jl")

export Planet_plane_hk, compute_ttv!
export Planet_plane
# export ttv_nplanet, ttv_wrapper, ttv_chisquare    

end
