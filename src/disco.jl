

function disco(X; method=:exact, n_sample=5000, alpha=1.0, rng=default_rng())
    if method == :exact
        return disco_exact(X, alpha)
    elseif method == :sample
        return disco_sample(X, n_sample, alpha, rng)
    else
        error("Unknown method='$(method)'")
    end
end

function disco_exact(X, alpha)

    within = 0.0
    for x in X
        within += pairdist_exact(x, alpha) * length(x) / 2
    end
    Xtot = vcat(X...)
    total = pairdist_exact(Xtot, alpha) * length(Xtot) / 2
    between = total - within

    return (total=total, between=between, within=within)
end

function disco_sample(X, n_sample, alpha, rng)

    within = 0.0
    for x in X
        within += pairdist_sample(x, n_sample, alpha, rng) * length(x) / 2
    end
    Xtot = vcat(X...)
    total = pairdist_sample(Xtot, n_sample, alpha, rng) * length(Xtot) / 2
    between = total - within

    return (total=total, between=between, within=within)
end
