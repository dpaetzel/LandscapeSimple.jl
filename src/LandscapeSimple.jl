module LandscapeSimple

using Random
using Sobol
using QuasiMonteCarlo

import Base: ∘

export mkscale_const,
    mkscale_discrete,
    mkscale_geo,
    mkscale_log,
    mkscale_minmax,
    mkscale_mix,
    configurations,
    sobols,
    sobols!

include("sobols.jl")

"""
    Rescaler(xmin, xmax, xminnew, xmaxnew)

A function that scales values from range [xmin, xmax] to new range [xminnew,
xmaxnew].

Essentially corresponds to

```
function rescale(
    x::T;
    xmin::T,
    xmax::T,
    xminnew::T,
    xmaxnew::T,
) where {T<:Number}
    return (x - xmin) / (xmax - xmin) * (xmaxnew - xminnew) + xminnew
end
```

but allows type stability.
"""
struct Rescaler{T<:Number}
    xmin::T
    xmax::T
    xminnew::T
    xmaxnew::T
    # TODO Consider to do input checking upon construction
end

"""
Always promotes the internals of the `Rescaler` and the incoming number to a
common type and transforms to the `Rescaler` type before output.
"""
function (rescaler::Rescaler{T}) where {T<:Number} end

function (rescaler::Rescaler{T})(x::AbstractArray) where {T<:Number}
    return rescaler.(x)
end

function (rescaler::Rescaler{T})(x) where {T<:Number}
    # TODO Consider to improve runtime performance by somehow not computing
    # derivatives wrt to parameters of this
    # TODO Consider to precompute this upon rescaler construction
    rescale =
        (rescaler.xmaxnew - rescaler.xminnew) / (rescaler.xmax - rescaler.xmin)
    out = (x - rescaler.xmin) * rescale + rescaler.xminnew
    return convert(T, out)
end

struct TypedScale{T}
    type::Type{T}
    transform
end

# We want to support composition.
function ∘(f, scale::TypedScale{T}) where {T<:Real}
    # Try to infer the output type of f applied to T.
    try
        sample_output = f(zero(T))
        OutputType = typeof(sample_output)
        return TypedScale(OutputType, f ∘ scale.transform)
    catch
        # We can always fallback to Any if we can't determine the type.
        return TypedScale(Any, f ∘ scale.transform)
    end
end

mkscale_const(x::T) where {T<:Real} = TypedScale(T, _ -> x)

"""
Generate a transformation that transforms a number in \$[0, 1]\$ to the interval
`[xmin, xmax]`.

Make sure that `T` is synchronized with `TSobol` of `configurations`.
"""
mkscale_minmax(xmin::T, xmax::T) where {T<:Real} =
    TypedScale(T, Rescaler(zero(T), one(T), xmin, xmax))

"""
Generate a transformation that transforms a number in \$[0, 1]\$ to the given set
of values, assigning them equal weights.
"""
function mkscale_discrete(vals::AbstractVector)
    function _scale_discrete(x)
        # Since we `ceil` below, we get index `0` if we put in `0` since
        # `ceil(0) == 0`.
        if x == 0
            return vals[1]
        else
            return vals[Int(ceil(x * length(vals)))]
        end
    end
    T = promote_type(typeof.(vals)...)
    return TypedScale(T, _scale_discrete)
end

"""
Generate a transformation that transforms a number in \$[0, 1]\$ to a geometric
progression between `xmin` and `xmax` with the given `base`.
"""
function mkscale_geo(xmin::T, xmax::T; base=10) where {T<:Real}
    if xmin < 0
        @warn "mkscale_geo was used with negative `xmin`, I hope you know " *
              "what you're doing …"
    end
    _scale_geo(x) =
        if base == 1
            xmin + x * (xmax - xmin)
        else
            xmin + (xmax - xmin) * (base^x - 1) / (base - 1)
        end
    # This should work but for now let's just require that types are uniform.
    # common_type = promote_type(typeof.((1, 1.0, 1.0f0))...)
    return TypedScale(T, _scale_geo)
end

"""
    mkscale_log(xmin, xmax; biaslow=1.0)

Construct a scaling function that maps a Sobol'/quasi-random draw `x ∈ [0,1]` to
a value in `[xmin, xmax]` on a log scale.

The mapping is log-uniform: small and large values are equally likely in
log-space.

The optional exponent `biaslow > 1` biases samples toward the lower end of the
range; `biaslow` < 1 biases toward the upper end; `biaslow = 1` is unbiased
log-uniform.
"""
function mkscale_log(xmin::T, xmax::T; biaslow=1.0) where {T<:Real}
    function _mkscale_log(x::T)
        return xmin * (xmax / xmin)^(x^biaslow)
    end

    return TypedScale(T, _mkscale_log)
end

"""
    mkscale_mix(scales; atol=1e-12, normalize=true)

Create a mixture scale from `(proportion => scale)` pairs.

- Proportions must be > 0. If `normalize=true` (default), they are normalised.
- Uses half-open bins [l,u) except the last bin [l,1].
- Returns a callable `_scale_mix(u::Real)` that maps `u ∈ [0,1)` to the mixed
  scale.

Each `scale` may be either:
  • a callable `scale(u)` expecting `u ∈ [0,1)`, or
  • an object with a `transform(u)` method/field.

Example:

```
julia> scales = [0.3 => mkscale_discrete([1,2,3]),
                 0.4 => mkscale_geo(0.1, 1.0),
                 0.3 => mkscale_minmax(0.2, 0.5)]
julia> mix = mkscale_mix(scales)
julia> y = mix(0.1234)
```
"""
function mkscale_mix(
    scales::Vector{P};
    atol::Real=1e-12,
    normalize::Bool=true,
) where {T,P<:Pair{<:Real,<:TypedScale{T}}}
    @assert !isempty(scales) "mkscale_mix requires at least one component"

    # Extract and validate proportions.
    raw_p = Float64[p for (p, _) in scales]
    @assert all(isfinite, raw_p) "mkscale_mix requires proportions to be finite"
    @assert all(p -> p > 0.0, raw_p) "mkscale_mix requires proportions to " *
                                     "be all positive"

    # Normalise (tolerantly) or enforce sum≈1.
    Σp = sum(raw_p)
    p = if normalize
        raw_p ./ Σp
    else
        @assert isapprox(Σp, 1.0; atol=atol) "mkscale_mix proportions " *
                                             "must sum to 1 " *
                                             "(|Σp-1| ≤ $atol)"
        raw_p
    end

    # Precompute cumulative edges and bin bounds.
    edges = cumsum(p)
    # Force exact 1.0.
    edges[end] = 1.0
    lbounds = [0.0; @view edges[1:(end - 1)]]
    ubounds = edges

    # Extract underlying scales.
    scs = [sc for (_, sc) in scales]

    # Pick bin: half-open [l,u) except last bin [l,1].
    @inline function pick_bin(u::Float64)
        # Clamp tiny overshoots to the first bin.
        if u ≤ 0.0 + eps(Float64)
            return 1
            # Clamp tiny overshoots to the last bin.
        elseif u ≥ 1.0 - eps(Float64)
            return length(lbounds)
        else
            # Last j with lbounds[j] ≤ u.
            #
            # “Return the index of the last value in v less than or equivalent
            # to x.  If x is less than all values in v the function returns
            # firstindex(v) - 1.”
            i = searchsortedlast(lbounds, u)
            # `searchsortedlast > length(edges)` if not found, clamp to last
            # element.
            return clamp(i, 1, length(lbounds))
        end
    end

    # Return the 1D mixer.
    function _scale_mix(u::Real)
        @assert isfinite(u) "u must be finite"
        # Clamp to [0, 1) to tolerate roundoff, but allow exactly 1.0 → last
        # bin.
        u64 = Float64(u)
        i = pick_bin(u64)
        l = lbounds[i]
        r = ubounds[i]
        w = r - l
        @assert w > 0.0 "internal error: zero-width bin"
        # Map u into local [0,1).
        v = (clamp(u64, l, r) - l) / w
        v = min(max(v, 0.0), 1.0 - eps(Float64))
        return scs[i].transform(v)
    end

    return TypedScale(T, _scale_mix)
end

# Note that since `sample_sobol` and `sample_sobol_scrambled` are deterministic,
# we can simply index into them and parallelize nicely at the Slurm level
# (however, we probably want to run a bunch of runs in each Julia process since
# Julia startup overhead is not 0).

"""
Convert a sample from \$[0, 1)^d\$ to a sample from the given hyperparameter
space.

# Arguments

- `space`: A hyperparameter space description. A `NamedTuple` mapping parameter
  names to `TypedScale`s.
- `v`: A sample from \$[0, 1)^d\$ (e.g. from a low discrepancy sequence).
"""
function _gethp(space, v)
    # Apply for each dimension in the space the transformation to the
    # corresponding value in `v`, then convert to the desired type of the
    # dimension. Annotate with the name of the dimension.
    return NamedTuple(
        name => convert(dim.type, map(dim.transform, v)) for
        (name, dim, v) in zip(keys(space), values(space), v)
    )
end

"""
Generate `2^m` configurations from the given space using the Sobol' sequence as
a source for quasi random numbers.

# Arguments

- `rng::AbstractRNG`: RNG used for scrambling the Sobol' sequence.
"""
function configurations(rng::AbstractRNG, space; m=9)
    # Have to filter for `<: Real` because `mkscale_discrete` may use
    # non-numeric target types.
    TsFloat = [v.type for v in values(space) if v.type <: AbstractFloat]
    T = if isempty(TsFloat)
        Float64
    else
        promote_type(TsFloat...)
    end
    sobols_ = sobols(rng, length(space), T; m=m)
    return configurations(sobols_, space)
end

configurations(space; m=9) = configurations(Random.Xoshiro(31), space; m=m)

"""
Having a sequence of low discrepancy vectors in ``[0, 1]^d`` (e.g. Sobol'
numbers), generate a configuration from the given d-dimensional space for each.
"""
function configurations(ldns, space)
    return _gethp.(Ref(space), eachcol(ldns))
end

end # module LandscapeSimple
