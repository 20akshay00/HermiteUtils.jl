module HermiteUtils

using FastGaussQuadrature

export hermite_polynomial, hermite_function
export hermite_polynomial_basis, hermite_function_basis
export hermite_product_integral

# ==============================================================================
# Core Recurrence
# ==============================================================================

"""
    _hermite_recurrence(f::F, n::Integer, x::Number, h0::Number) where {F}

Internal function that computes the normalized Hermite recurrence up to degree `n`.
Applies the callback function `f(k, h)` at each degree `k` with value `h`. 
"""
function _hermite_recurrence(f::F, n::Integer, x::Number, h0::Number) where {F}
    n < 0 && throw(DomainError(n, "Degree n must be non-negative"))

    f(0, h0)
    n == 0 && return

    h1 = sqrt(2) * x * h0
    f(1, h1)
    n == 1 && return

    hm2, hm1 = h0, h1
    for k in 2:n
        h = sqrt(2 / k) * x * hm1 - sqrt((k - 1) / k) * hm2
        f(k, h)
        hm2, hm1 = hm1, h
    end
end

# core wrappers to extract the final value or fill an array
function _evaluate_hermite(n::Integer, x::Number, h0::Number)
    out = Ref{typeof(h0 * x)}()
    _hermite_recurrence(n, x, h0) do k, h
        out[] = h
    end
    return out[]
end

function _evaluate_hermite_basis!(out::AbstractVector, n::Integer, x::Number, h0::Number)
    _hermite_recurrence(n, x, h0) do k, h
        out[k+1] = h
    end
    return out
end

# ==============================================================================
# Normalized Hermite Polynomials
# ==============================================================================

"""
    hermite_polynomial(n::Integer, x::Number)
    hermite_polynomial(n::Integer)

Evaluate the normalized physicist's Hermite polynomial of degree `n` at `x`.
If `x` is omitted, returns a function `x -> hermite_polynomial(n, x)`.
"""
function hermite_polynomial(n::Integer, x::Number)
    T = typeof(float(x))
    h0 = T(π)^(-0.25)
    return _evaluate_hermite(n, x, h0)
end

hermite_polynomial(n::Integer) = x -> hermite_polynomial(n, x)

"""
    hermite_polynomial_basis(n::Integer, x::Number)

Return a vector of normalized Hermite polynomials from degree 0 up to `n` evaluated at `x`.
"""
function hermite_polynomial_basis(n::Integer, x::Number)
    T = typeof(float(x))
    out = Vector{T}(undef, n + 1)
    h0 = T(π)^(-0.25)
    return _evaluate_hermite_basis!(out, n, x, h0)
end


# ==============================================================================
# Normalized Hermite Functions
# ==============================================================================

"""
    hermite_function(n::Integer, x::Number)
    hermite_function(n::Integer)

Evaluate the normalized Hermite function (includes the `exp(-x^2 / 2)` term) of 
degree `n` at `x`. If `x` is omitted, returns a function `x -> hermite_function(n, x)`.
"""
function hermite_function(n::Integer, x::Number)
    T = typeof(float(x))
    h0 = T(π)^(-0.25) * exp(-x^2 / 2)
    return _evaluate_hermite(n, x, h0)
end

hermite_function(n::Integer) = x -> hermite_function(n, x)

"""
    hermite_function_basis(n::Integer, x::Number)

Return a vector of normalized Hermite functions from degree 0 up to `n` evaluated at `x`.
"""
function hermite_function_basis(n::Integer, x::Number)
    T = typeof(float(x))
    out = Vector{T}(undef, n + 1)
    h0 = T(π)^(-0.25) * exp(-x^2 / 2)
    return _evaluate_hermite_basis!(out, n, x, h0)
end


# ==============================================================================
# Integrals
# ==============================================================================

"""
    hermite_product_integral(ns::Tuple{Vararg{Integer}}; quad_rule=nothing)

Analytically compute the integral of the product of `m` Hermite functions with 
degrees specified by the tuple `ns`. 

Uses Gauss-Hermite quadrature, which is exact for this integration.
If doing multiple integrals in a tight loop, you can precompute and pass 
`quad_rule = gausshermite(N)` to avoid recomputing quadrature nodes.

# Example
```julia
# Integral of ψ_0(x) * ψ_1(x) * ψ_2(x)
hermite_product_integral((0, 1, 2))
"""

function hermite_product_integral(ns::Tuple{Vararg{Integer}}; quad_rule=nothing)
    isodd(sum(ns)) && return 0.0
    nmax = maximum(ns)
    m = length(ns)

    if quad_rule === nothing
        N = ceil(Int, (sum(ns) + 1) / 2)
        quad_rule = gausshermite(N)
    end

    xs, ws = quad_rule
    s = sqrt(2 / m)

    T = typeof(float(s * xs[1]))
    basis = Vector{T}(undef, nmax + 1)
    h0 = T(π)^(-0.25) 

    integral = 0.0
    for i in eachindex(ws)
        x_eval = s * xs[i]

        _evaluate_hermite_basis!(basis, nmax, x_eval, h0)

        prod_val = 1.0
        for n in ns
            prod_val *= basis[n+1]
        end

        integral += ws[i] * prod_val
    end

    return s * integral
end

end