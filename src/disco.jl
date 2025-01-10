

function disco(X; method=:exact, n_sample=5000, alpha=1.0)
    if method == :exact
        return disco_exact(X, alpha)
    elseif method == :sample
        return disco_sample(X, n_sample, alpha)
    else
        error("Unknown method='$(method)'")
    end
end

function disco_exact(X, alpha)

    within = 0.0
    for x in X
        within += pairdist_exact(x, x, alpha) * length(x) / 2
    end
    Xtot = vcat(X...)
    total = pairdist_exact(Xtot, Xtot, alpha) * length(Xtot) / 2
    between = total - within

    return (total=total, between=between, within=within)
end

function disco_sample(X, n_sample, alpha)

    within = 0.0
    for x in X
        within += pairdist_sample(x, x, n_sample, alpha) * length(x) / 2
    end
    Xtot = vcat(X...)
    total = pairdist_sample(Xtot, Xtot, n_sample, alpha) * length(Xtot) / 2
    between = total - within

    return (total=total, between=between, within=within)
end
