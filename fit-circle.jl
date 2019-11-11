module Circle

using LinearAlgebra

function generate_test(filename, count, origin, radius, noise)
    open(filename, "w") do f
        for _ in 1:count
            alpha = rand() * 2pi
            p = origin + [cos(alpha), sin(alpha)] * radius + (rand(2) .- 0.5) * noise
            println(f, "$(p[1]) $(p[2])")
        end
    end
end

function generate_test_3d(filename, count, origin, radius, noise)
    open(filename, "w") do f
        for _ in 1:count
            theta = rand() * 2pi
            phi = rand() * 2pi
            p = origin + [sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)] * radius +
                (rand(3) .- 0.5) * noise
            println(f, "$(p[1]) $(p[2]) $(p[3])")
        end
    end
end

function read_data(filename)
    data = []
    open(filename) do f
        for line in eachline(f)
            push!(data, map(s -> parse(Float64, s), split(line)))
        end
    end
    data
end

function write_result(filename, count, origin, radius)
    open(filename, "w") do f
        for i in 1:count
            alpha = (i - 1) / (count - 1) * 2pi
            p = origin + [cos(alpha), sin(alpha)] * radius
            println(f, "$(p[1]) $(p[2])")
        end
    end
end

function intersect_circles(p1, p2, r)
    d = p2 - p1
    n = [d[2], -d[1]]
    x = 0.5
    y2 = (r / norm(d))^2 - 0.25
    (y2 < 0 || isnan(y2) || isinf(y2)) && return nothing, nothing
    y = sqrt(y2)
    p1 + d * x + n * y, p1 + d * x - n * y
end

function guess_circle(p1, p2, p3, r)
    a1, a2 = intersect_circles(p1, p2, r)
    b1, b2 = intersect_circles(p1, p3, r)
    any(x -> x === nothing, [a1, a2, b1, b2]) && return nothing
    local best
    best_dist = Inf
    for (a, b) in [(a1, b1), (a1, b2), (a2, b1), (a2, b2)]
        d = norm(a - b)
        if d < best_dist
            best = (a + b) / 2
            best_dist = d
        end
    end
    best
end

function fit_circle(data, r, tolerance, max_iteration)
    result = nothing
    while result === nothing
        result = guess_circle(rand(data), rand(data), rand(data), r)
    end
    for i in 1:max_iteration
        guess = nothing
        while guess === nothing
            guess = guess_circle(rand(data), rand(data), rand(data), r)
        end
        next = (result * i + guess) / (i + 1)
        norm(next - result) < tolerance && return next
        result = next
    end
    @warn "Exited with maximum iteration"
    result
end

function lsq_fit(X, R, tolerance = 1.0e-7, max_iteration = 1000, λ = 0.1)
    n = length(X[1])            # dimension
    m = length(X)
    a = sum(X) / m
    J = Array{Float64}(undef, n * m, n)
    rhs = Array{Float64}(undef, n * m)
    id = Array{Float64}(I, n, n)
    for _ in 1:max_iteration
        for i in 1:m
            v = X[i] - a
            d = norm(v)
            J[(i-1)*n+1:i*n,:] = id - (id - v * transpose(v) / d^2) * R / d
            rhs[(i-1)*n+1:i*n] = v / d * (d - R)
        end
        Δa = J \ rhs
        a += λ * Δa
        norm(Δa) < tolerance && return a
    end
    @warn "Exited with maximum iteration"
    a
end

function lsq_fit_noradius(X, tolerance = 1.0e-7, max_iteration = 1000, λ = 0.1)
    n = length(X[1]) + 1
    m = length(X)
    Xc = sum(X) / m
    R = sqrt(sum(Xi -> norm(Xi - Xc)^2, X) / m)
    a = vcat(R, Xc)
    J = Array{Float64}(undef, n * m, n)
    rhs = Array{Float64}(undef, n * m)
    dRda = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
    dXda = [0 1 0 0; 0 0 1 0; 0 0 1 0; 0 0 0 1]
    for _ in 1:max_iteration
        R = a[1]
        Xc = a[2:n]
        for i in 1:m
            v = X[i] - Xc
            d = norm(v)
            J[(i-1)*n+1:i*n,:] = dXda + v / d * dRda - R / d * (I - dot(v, v) / d^2) * dXda
            rhs[(i-1)*n+1:i*n] = v / d * (d - R)
        end
        Δa = J \ rhs
        a += λ * Δa
        norm(Δa) < tolerance && return a
    end
    @warn "Exited with maximum iteration"
    a[1], a[2:n]
end

function test()
    r = 10.0
    origin = [2, 5.]
    noise = 1.0
    generate_test("/tmp/points", 100, origin, r, noise)
    write_result("/tmp/original", 100, origin, r)
    data = read_data("/tmp/points")
    @time p = fit_circle(data, r, 1.0e-7, 1000)
    println("[Guess] Error: $(norm(p - origin))")
    @time p = lsq_fit(data, r, 1.0e-7, 1000)
    write_result("/tmp/fit", 100, p, r)
    println("[ LSQ ] Error: $(norm(p - origin))")
    @time R, p = lsq_fit_noradius(data, 1.0e-7, 1000)
    write_result("/tmp/fit", 100, p, R)
    println("[LSQwo] Error: $(norm(p - origin)) and $(abs(R - r))")
end

function test3d()
    r = 10.0
    origin = [2, 5, 3.]
    noise = 2.0
    generate_test_3d("/tmp/points", 1000, origin, r, noise)
    data = read_data("/tmp/points")
    @time p = lsq_fit(data, r, 1.0e-7, 1000)
    println("Error: $(norm(p - origin))")
end

end
