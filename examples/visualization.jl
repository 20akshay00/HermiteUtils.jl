using HermiteUtils
using Plots
using FastGaussQuadrature

# plotting polynomials and functions
# ==============================================================================

begin
    theme(:dao)
    x_poly = range(-2.5, 2.5, length=300)
    x_func = range(-6.0, 6.0, length=300)
    n_max = 10
    grad = cgrad(:rainbow)
end;

begin
    # plot polynomials shifted vertically for clarity
    p1 = plot(title="Orthonormal Hermite polynomials", legend=:right)
    for n in 0:n_max
        # compute basis at each x to build the curve
        y = [hermite_polynomial(n, xi) for xi in x_poly]
        # add a vertical offset of 3.0 per degree to stack them
        plot!(p1, x_poly, y .+ 2n, label="", lw=1.5, c=grad[n/n_max])
    end
    savefig("./examples/hermite_polynomials.png")
    p1
end

begin
    # plot hermite functions shifted vertically
    p2 = plot(title="Orthonormal Hermite functions", legend=:right)
    for n in 0:n_max
        # using the basis function is more efficient for full evaluation
        y = [hermite_function(n, xi) for xi in x_func]
        # add a vertical offset of 1.5 per degree
        plot!(p2, x_func, y .+ 1.5n, label="", lw=1.5, c=grad[n/n_max])
    end
    savefig("./examples/hermite_functions.png")
    p2
end

# four-integral decay
# ==============================================================================
begin
    n_limit = 4 # indices 0:4
    results = []

    # precompute quadrature rule for efficiency
    # max degree in the product is sum(ns) = 4 * n_limit
    quad = gausshermite(ceil(Int, (4 * n_limit + 1) / 2))

    # iterate through unique combinations only
    for i in 0:n_limit
        for j in i:n_limit
            for k in j:n_limit
                for l in k:n_limit
                    val = hermite_product_integral((i, j, k, l); quad_rule=quad)
                    push!(results, ("($i,$j,$k,$l)", abs(val)))
                end
            end
        end
    end

    # sort results by magnitude to show the decay of coefficients
    sort!(results, by=x -> x[2], rev=true)
end

begin
    nvals = 34
    labels = [r[1] for r in results][1:nvals]
    values = [r[2] for r in results][1:nvals]

    p3 = bar(values,
        xticks=(1:length(labels), labels),
        xrotation=90,
        title="\nDecay of 4-Hermite integrals",
        ylabel="Magnitude",
        xlabel="Indices (i,j,k,l)",
        size=(1000, 500),
        tickfont=font(7),
        ylims=[-0.01, 0.45],
        xlims=[0, nvals],
        legend=nothing,
        margin=10Plots.mm
    )
    savefig("./examples/hermite_four_integrals.png")
    display(p3)
end