module LandscapeSimple

using Random
using Sobol
using QuasiMonteCarlo

export mkscale_discrete,
    mkscale_geo, mkscale_minmax, configurations, gethp, sobols, sobols!

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

function (rescaler::Rescaler{T})(x::AbstractArray{T}) where {T<:Number}
    return rescaler.(x)
end

function (rescaler::Rescaler{T})(x::T) where {T<:Number}
    # TODO Consider to improve runtime performance by somehow not computing
    # derivatives wrt to parameters of this
    # TODO Consider to precompute this upon rescaler construction
    rescale =
        (rescaler.xmaxnew - rescaler.xminnew) / (rescaler.xmax - rescaler.xmin)
    return (x - rescaler.xmin) * rescale + rescaler.xminnew
end

"""
Generate a transformation that transforms a number in \$[0, 1]\$ to the interval
`[xmin, xmax]`.
"""
mkscale_minmax(xmin, xmax) = Rescaler(0.0f0, 1.0f0, xmin, xmax)

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
    return _scale_discrete
end

"""
Generate a transformation that transforms a number in \$[0, 1]\$ to a geometric
progression between `xmin` and `xmax` with the given `base`.
"""
function mkscale_geo(xmin, xmax; base=10)
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
    return _scale_geo
end

# Note that since `sample_sobol` and `sample_sobol_scrambled` are deterministic,
# we can simply index into them and parallelize nicely at the Slurm level
# (however, we probably want to run a bunch of runs in each Julia process since
# Julia startup overhead is not 0).

"""
Convert a sample from \$[0, 1)^d\$ to a sample from the given hyperparameter
space.
"""
function gethp(space, v)
    return (; zip(keys(space), map.(values(space), v))...)
end

"""
Generate `2^m` configurations from the given space using the Sobol' sequence as
a source for quasi random numbers.

# Arguments

- `rng::AbstractRNG`: RNG used for scrambling the Sobol' sequence.
"""
function configurations(rng::AbstractRNG, space; m=9)
    sobols_ = sobols(rng, length(space); m=m)
    return configurations(sobols_, space)
end

configurations(space; m=9) = configurations(Random.Xoshiro(31), space; m=m)

"""
Having a sequence of low discrepancy vectors in ``[0, 1]^d`` (e.g. Sobol'
numbers), generate a configuration from the given d-dimensional space for each.
"""
function configurations(ldns, space)
    return gethp.(Ref(space), eachcol(ldns))
end

end # module LandscapeSimple
