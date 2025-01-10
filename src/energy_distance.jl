
function energy_distance(X, Y; method=:exact, n_sample=5000, alpha=1.0)
    if method == :exact
        return energy_distance_exact(X, Y, alpha)
    elseif method == :sample
        return energy_distance_sample(X, Y, n_sample, alpha)
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

function pairdist_sample(X, Y, n_sample, alpha)

    # Better to use exact calculation in this case
    if n_sample >= length(X) * length(Y)
        return pairdist_exact(X, Y, alpha)
    end

    d = 0.0
    for _ in 1:n_sample
        d += norm(sample(X) - sample(Y))^alpha
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

function energy_distance_sample(X, Y, alpha, n_sample)
    dxy = pairdist_sample(X, Y, alpha, n_sample)
    dxx = pairdist_sample(X, X, alpha, n_sample)
    dyy = pairdist_sample(Y, Y, alpha, n_sample)
    return 2*dxy - dxx - dyy
end
