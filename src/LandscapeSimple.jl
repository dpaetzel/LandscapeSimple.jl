module LandscapeSimple

using Random
using QuasiMonteCarlo

export mkscale_discrete, mkscale_minmax, configurations, gethp, sobols

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
Get \$2^m\$ scrambled Sobol' numbers for the given dimensionality.

Note that we actually generate \$2^{m+1}\$ numbers using `Sobol.jl` and then
skip the first \$2^m\$ to account for the fact that the numbers are otherwise
not balanced due to `Sobol.jl` not starting sequences with \$(0)^d\$.

# Arguments

- `rng::AbstractRNG`: RNG used for scrambling the Sobol' sequence.

# Example

```
julia> UnicodePlots.scatterplot(eachrow(sobols(2))...)
     ┌────────────────────────────────────────┐
   1 │⡀⠀⠂⠁⠄⡐⢀⠠⠀⠈⠂⠀⡀⠄⠁⢀⠐⠈⠀⠠⠠⠈⠐⠀⠠⠁⠀⠄⡀⢁⠈⠠⢀⠀⠐⠂⠀⠁⠂⡀│
     │⠐⠈⠐⠀⠠⠀⠀⠄⣀⠁⠈⠠⢀⠀⠈⠂⠀⠁⠂⠄⠄⠀⠂⠁⠂⠐⣀⠠⠀⠀⠁⠀⡀⠄⠁⠠⡐⠈⠀⠠│
     │⡀⠢⠀⠂⠈⠄⢀⠁⠠⠀⠀⣁⠀⠀⠠⠐⠠⠀⠊⠀⠈⠐⠀⠌⢀⠀⠄⠈⡀⠠⠂⠀⠀⣈⠀⠀⠐⠀⠄⠐│
     │⠀⠀⡂⠨⠀⢀⠂⠀⠁⠠⡂⠀⠀⠡⠀⠀⢃⠠⠀⠈⠀⠆⢐⠀⠀⠠⠈⠀⠐⡀⠀⠌⠀⠀⢂⠁⠀⠄⡘⠀│
     │⠁⠀⢀⠁⠄⡀⠐⠠⠀⡈⠀⠀⠂⠄⠡⠐⠀⠈⡀⠐⠠⠈⡀⠀⠐⢁⠀⠄⠂⢀⠌⠠⠐⠀⡀⠄⢀⠁⠀⠂│
     │⠐⡈⠀⠀⠐⠀⠠⠄⡂⠀⠈⠠⠐⡀⠀⡂⠀⠁⢀⠁⠁⠀⠀⢁⠂⠀⢐⠠⠄⠀⠀⢀⠂⠤⠀⠈⠀⠈⠀⢐│
     │⡀⠁⠐⢀⠈⠁⡀⠠⠠⠀⠀⠂⡀⠠⠂⠀⠀⠁⠐⠠⠁⡀⠂⠈⡀⠀⠄⠄⢀⠂⠐⠄⢀⠐⠀⠠⠃⠈⢀⠀│
     │⠀⠐⠄⡁⠀⠠⠂⠀⠀⢐⠆⠀⠀⠈⡀⢀⠱⢀⠁⠀⠀⢁⠠⠂⠀⡁⠀⠀⠐⠄⢀⠁⠀⠀⠰⠀⠈⡀⠔⠀│
     │⠁⢀⠠⠀⢁⡀⠈⠐⠀⠂⠠⠀⠁⠂⠠⠀⠀⡀⠄⠐⡈⠀⠄⡀⠈⠐⠀⠂⠁⠠⠄⠐⠈⠀⠄⢂⠠⢀⠀⠁│
     │⠐⡄⠀⠀⠠⠀⠀⠂⠡⠀⢀⠒⠈⠄⠀⠅⠀⡀⢀⠁⡀⠀⠀⢠⠁⠀⠌⠐⠀⠀⠀⠀⠅⠒⠀⠀⡀⢀⠀⠨│
     │⢀⠀⠈⢐⠀⠈⠄⠐⠈⢀⠀⠀⠄⢐⠂⠠⠀⢀⠈⠄⠐⡀⠁⠀⠄⡈⠀⠂⠠⠁⠐⡂⠠⠀⠁⠠⠁⡀⠀⠄│
     │⠀⠀⠅⠂⡀⠠⠁⢀⠀⠐⠅⠀⠀⡀⠂⢀⠡⠐⠀⢀⠠⠐⠠⠁⠀⠂⠀⡀⠈⠄⠐⢀⠀⠀⠨⡀⠀⠂⠌⡀│
     │⠡⠐⠠⠀⢀⠂⠈⡀⠀⠂⠈⢀⠁⠀⠠⠄⠀⠂⠄⡀⡀⠀⠄⠂⡈⠐⠀⢀⠁⠐⠄⠀⠈⡀⠁⢀⠠⠐⠀⠠│
     │⠠⡀⠀⠄⠂⠀⠀⠂⠩⠀⠀⡂⠈⠀⢀⠠⢀⠀⠠⠂⠈⠀⠄⢠⠀⠀⠍⠐⠀⠀⡀⠀⠁⢐⠀⠐⠄⠀⡀⠂│
   0 │⠂⠀⠄⢐⠀⠡⠁⠐⠀⡀⠠⠈⠀⠐⡀⠀⠅⢀⠀⠐⠀⡂⠠⠀⠐⢀⠀⠂⠈⡈⢀⠂⠀⠁⠄⠄⠀⡀⠨⠀│
     └────────────────────────────────────────┘
     ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀1⠀
```
"""
function sobols(rng, dims; m=9)
    # Sobol sequence only has nice properties if generating `2^m` points for
    # some `m`.
    N = 2^m
    # This is fixed for the Sobol' sequence.
    b = 2

    # Draw `2N` and throw away the first `N` to work around `Sobol.jl` not
    # including `(0)^d`; see
    # https://github.com/JuliaMath/Sobol.jl/issues/31#issue-1328916662
    # (QuasiMonteCarlo.jl is based on Sobol.jl).
    sample_sobol =
        QuasiMonteCarlo.sample(2 * N, dims, SobolSample(; R=NoRand()), Float32)
    # Throw away first `N` points (see comment above on `(0)^dims` not being
    # included).
    sample_sobol = sample_sobol[:, (N + 1):end]

    # TODO Consider to additionally do Digital Shift Scrambling like SciPy.
    # `pad=32` is the default as of 2024-06-26.
    # `MatousekScramble` is likely what SciPy means by linear matrix scramble (LMS).
    sample_sobol_scrambled =
        randomize(sample_sobol, MatousekScramble(; base=b, pad=32, rng=rng))

    return sample_sobol_scrambled
end

sobols(dims; m=9) = sobols(Random.Xoshiro(31), dims; m=m)

"""
Generate `2^m` configurations from the given space using the Sobol' sequence as
a source for quasi random numbers.

# Arguments

- `rng::AbstractRNG`: RNG used for scrambling the Sobol' sequence.
"""
function configurations(rng::AbstractRNG, space; m=9)
    return gethp.(Ref(space), eachcol(sobols(rng, length(space); m=m)))
end

configurations(space; m=9) = configurations(Random.Xoshiro(31), space; m=m)

end # module LandscapeSimple
