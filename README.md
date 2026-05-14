# HermiteUtils.jl

A simple utility package for working with normalized physicist's Hermite polynomials and Hermite functions.

## Installation

Since this is a personal package, you can add it via its local path or its GitHub URL:

```julia
using Pkg
Pkg.add("https://github.com/20akshay00/HermiteUtils.jl")
```

## Usage

The core of the package is a single recurrence relation:
$h_{k} = \sqrt{\frac{2}{k}} x h_{k-1} - \sqrt{\frac{k-1}{k}} h_{k-2}$
Building upon this, it provides a handful of wrappers that are easier to use.

### Hermite Polynomials
These functions evaluate the normalized physicist's Hermite polynomials $\tilde{H}_n(x)$. Note that this includes the normalization factors such that $\int_{-\infty}^{\infty} e^{-x^2} \tilde{H}_m(x) \tilde{H}_n(x) dx = \delta_{mn}$.

```julia
# evaluate a single polynomial of degree n at x
val = hermite_polynomial(3, 1.5)

# return a closure for a specific degree
H3 = hermite_polynomial(3)
y = H3(1.5)

# evaluate all polynomials from degree 0 to n at x
# this returns a vector of length n + 1
basis = hermite_polynomial_basis(5, 0.5)
```

### Hermite Functions
These functions evaluate the normalized Hermite functions $\psi_n(x) = \tilde{H}_n(x) e^{-x^2/2}$.

```julia
# evaluate a single hermite function at x
val = hermite_function(2, 0.5)

# get a function closure for a specific degree
ψ2 = hermite_function(2)
y = ψ2(0.5)

# evaluate the basis from degree 0 to n
Φ = hermite_function_basis(10, -1.2)
```

### Product Integrals
The package can also compute the integral of an arbitrary product of Hermite functions: $\int_{-\infty}^{\infty} \psi_{n_1}(x) \psi_{n_2}(x) \dots \psi_{n_m}(x) dx$. This is performed using Gauss-Hermite quadrature, which is analytically exact for these products.

```julia
# verifying normalization
hermite_product_integral((2, 2)) ≈ 1.0

# compute the integral of a four-function product
res = hermite_product_integral((0, 1, 1, 2))

# use a precomputed quadrature rule for high-performance loops
using FastGaussQuadrature
quad = gausshermite(20)
for i in 1:10
    hermite_product_integral((i, i, 0, 0); quad_rule=quad)
end
```
