module EnergyStats

    using StatsBase, LinearAlgebra
    import Random: default_rng
    export energy_distance, disco, dcov, dcor

    include("energy_distance.jl")
    include("disco.jl")
    include("dcov.jl")

end
