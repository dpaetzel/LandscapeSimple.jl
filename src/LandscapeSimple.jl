module LandscapeSimple

using Random
using Sobol
using QuasiMonteCarlo

import Base: ∘

export mkscale_discrete,
    mkscale_geo, mkscale_minmax, configurations, sobols, sobols!

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
