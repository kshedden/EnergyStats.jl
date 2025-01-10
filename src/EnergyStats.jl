module EnergyStats

    using StatsBase, LinearAlgebra

    export energy_distance, disco, dcov, dcor

    include("energy_distance.jl")
    include("disco.jl")
    include("dcov.jl")

end
