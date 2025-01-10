using EnergyStats
using Test
using StableRNGs

@testset "Energy distance" begin

    n, p = 100, 5
    rng = StableRNG(123)

    for i in 1:10
        X = [randn(rng, p) for _ in 1:n]
        m = (i - 1) / 10
        Y = [m .+ randn(rng, p) for _ in 1:2*n]

        d1 = energy_distance(X, Y; method=:exact)
        d2 = energy_distance(X, Y; method=:sample, n_sample=5000)
        @test isapprox(d1, d2, rtol=0.1, atol=0.1)
    end
end

@testset "Disco" begin

    n, m, p = 20, 20, 5
    rng = StableRNG(123)

    for i in 1:10
        X = [[randn(rng, p) for _ in 1:n] for _ in 1:m]
        d1 = disco(X; method=:exact)
        d2 = disco(X; method=:sample)
        @test isapprox(d1.total, d2.total, rtol=0.1, atol=0.1)
        @test isapprox(d1.within, d2.within, rtol=0.1, atol=0.1)
        @test isapprox(d1.between, d2.between, rtol=0.3, atol=0.1)
    end

end

@testset "Distance covariance and correlation" begin

    n, p = 20, 5
    rng = StableRNG(123)

    # Test the null case
    for debias in [false, true]
        for i in 1:10
            X = [randn(rng, p) for _ in 1:n]
            Y = [randn(rng, p) for _ in 1:n]
            d1 = dcov(X, Y; debias=debias)
            println(d1)
        end
    end

end
