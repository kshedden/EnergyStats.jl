
function energy_distance(X, Y; method=:exact, n_sample=5000, alpha=1.0, rng=default_rng())
    if method == :exact
        return energy_distance_exact(X, Y, alpha)
    elseif method == :sample
        return energy_distance_sample(X, Y, n_sample, alpha, rng)
    else
        error("Unknown method='$(method)'")
    end
end

function pairdist_exact(X, Y, alpha)
    d = 0.0
    for x in X
        for y in Y
            d += norm(x - y)^alpha
        end
    end
    m = length(X) * length(Y)
    d /= m
    return d
end

function pairdist_exact(X, alpha)
    d = 0.0
    n = length(X)
    for i in 1:n
        for j in 1:(i-1)
            d += norm(X[i] - X[j])^alpha
        end
    end
    m = Int(n*(n-1)/2)
    d /= m
    return d
end

function pairdist_sample(X, Y, n_sample, alpha, rng)

    # Better to use exact calculation in this case
    if n_sample >= length(X) * length(Y)
        return pairdist_exact(X, Y, alpha)
    end

    d = 0.0
    for _ in 1:n_sample
        d += norm(sample(rng, X) - sample(rng, Y))^alpha
    end
    d /= n_sample
    return d
end

function pairdist_sample(X, n_sample, alpha, rng)

    n = length(X)
    m = Int(n*(n-1)/2)

    # Better to use exact calculation in this case
    if n_sample >= m
        return pairdist_exact(X, alpha)
    end

    d = 0.0
    for _ in 1:n_sample
        ii = sample(rng, 1:n, 2; replace=false)
        d += norm(X[ii[1]] - X[ii[2]])^alpha
    end
    d /= n_sample
    return d
end

function energy_distance_exact(X, Y, alpha)
    dxy = pairdist_exact(X, Y, alpha)
    dxx = pairdist_exact(X, X, alpha)
    dyy = pairdist_exact(Y, Y, alpha)
    return 2*dxy - dxx - dyy
end

function energy_distance_sample(X, Y, alpha, n_sample, rng)
    dxy = pairdist_sample(X, Y, alpha, n_sample, rng)
    dxx = pairdist_sample(X, X, alpha, n_sample, rng)
    dyy = pairdist_sample(Y, Y, alpha, n_sample, rng)
    return 2*dxy - dxx - dyy
end
