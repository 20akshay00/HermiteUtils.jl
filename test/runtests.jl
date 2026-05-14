using Test
using HermiteUtils
using FastGaussQuadrature

using SpecialFunctions, HypergeometricFunctions, StaticArrays

# analytic baseline
function hermite_four_integral(i, j, k, l; max_level=typemax(Int))
    all(0 .≤ (i, j, k, l) .≤ max_level) || return 0.0
    iseven(i + j + k + l) || return 0.0

    # enforce correct symmetry and k ≥ l
    i, j, k, l = sort(SVector(i, j, k, l); rev=true)

    a = i + j - k + l + 1
    b = i - j + k - l + 1
    c = -i + j + k - l + 1
    d = -i - j + k - l + 1

    p = sqrt(2 * gamma(i + 1) * gamma(j + 1) * gamma(k + 1) * gamma(l + 1)) * pi^2
    q = gamma(a / 2) * gamma(b / 2) * gamma(c / 2) * gamma(k + 1) / gamma(k - l + 1)

    f1 = _₃F₂(float(-l), b / 2, c / 2, float(1 + k - l), d / 2, 1.0)

    if isnan(f1)
        # workaround for some issues with large arguments
        fp = _₃F₂(float(-l), b / 2, c / 2, float(1 + k - l), d / 2, 1.0 + eps())
        fm = _₃F₂(float(-l), b / 2, c / 2, float(1 + k - l), d / 2, 1.0 - eps())
        f = (fp + fm) / 2
    else
        f = f1
    end

    return f * q / p
end

@testset "HermiteUtils.jl" begin

    @testset "Basic checks" begin
        @test hermite_function(0, 0.0) ≈ π^(-0.25)
        @test hermite_polynomial(0, 0.0) ≈ π^(-0.25)

        for n in (1, 3, 5)
            @test hermite_function(n, 0.0) == 0.0
            @test hermite_polynomial(n, 0.0) == 0.0
        end

        ψ2 = hermite_function(2)
        @test ψ2(1.5) == hermite_function(2, 1.5)
    end

    @testset "Basis vs. scalar consistency" begin
        x = 1.23
        nmax = 5

        # functions
        f_basis = hermite_function_basis(nmax, x)
        @test length(f_basis) == nmax + 1
        for k in 0:nmax
            @test f_basis[k+1] ≈ hermite_function(k, x)
        end

        # polynomials
        p_basis = hermite_polynomial_basis(nmax, x)
        @test length(p_basis) == nmax + 1
        for k in 0:nmax
            @test p_basis[k+1] ≈ hermite_polynomial(k, x)
        end
    end

    @testset "Orthogonality (2-product integrals)" begin
        # ∫ ψ_m(x) ψ_n(x) dx = δ_mn
        for m in 0:4
            for n in 0:4
                ans = hermite_product_integral((m, n))
                if m == n
                    @test ans ≈ 1.0 atol = 1e-14
                else
                    @test ans ≈ 0.0 atol = 1e-14
                end
            end
        end
    end

    @testset "Higher-order integrals & symmetries" begin
        # odd sums should be exactly 0.0
        @test hermite_product_integral((1, 2, 2)) == 0.0
        @test hermite_product_integral((0, 1, 3, 5)) == 0.0

        # permutation symmetry
        val1 = hermite_product_integral((1, 2, 3))
        val2 = hermite_product_integral((3, 1, 2))
        @test val1 ≈ val2 atol = 1e-14

        # precomputed quadrature rule
        quad = gausshermite(10)
        val_default = hermite_product_integral((2, 2, 2, 2))
        val_precomp = hermite_product_integral((2, 2, 2, 2); quad_rule=quad)
        @test val_default ≈ val_precomp atol = 1e-14
    end

    @testset "Analytic 4-product check" begin
        for i in 0:4, j in 0:4, k in 0:4, l in 0:4
            val_quad = hermite_product_integral((i, j, k, l))
            val_analytic = hermite_four_integral(i, j, k, l)
            @test val_quad ≈ val_analytic atol = 1e-10
        end
    end
end