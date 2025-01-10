
function double_center(X)
    D = [norm(x1 - x2) for x1 in X, x2 in X]
    Dc = mean(D; dims=1)
    Dr = mean(D; dims=2)
    Dm = mean(D)
    D .-= Dr
    D .-= Dc
    D .+= Dm
    return D
end

function double_center_bc(X)
    n = length(X)
    D = [norm(x1 - x2) for x1 in X, x2 in X]
    Dc = sum(D; dims=1)
    Dr = sum(D; dims=2)
    Dm = sum(Dc)
    D .-= Dr ./ (n - 2)
    D .-= Dc ./ (n - 2)
    D .+= Dm ./ ((n - 1) * (n - 2))
    return D
end

function odsum(A, B)
    s = 0.0
    n = size(A, 1)
    @assert all(size(A) .== size(B))
    for i in 1:n
        for j in 1:n
            if i != j
                s += A[i, j] * B[i, j]
            end
        end
    end
    return s
end

function dcov(X, Y; debias=false)

    if length(X) != length(Y)
        @error("Lengths of X and Y must be equal")
    end

    n = length(X)

    if debias
        Dx = double_center_bc(X)
        Dy = double_center_bc(Y)
        return odsum(Dx, Dy) / (n * (n - 3))
    else
        Dx = double_center(X)
        Dy = double_center(Y)
        return sqrt(sum(Dx .* Dy) / n^2)
    end

end


function dcor(X, Y; debias=false)

    dxy = dcov(X, Y; debias=debias)
    dxx = dcov(X, X; debias=debias)
    dyy = dcov(Y, Y; debias=debias)

    return dxy / sqrt(dxx * dyy)
end
