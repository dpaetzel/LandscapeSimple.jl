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
    if N <= 0
        throw(
            ArgumentError(
                "The `m` you chose yielded N = 2^m <= 0 probably " *
                "because it is too large",
            ),
        )
    end

    # This is the default for the Matoušek scramble.
    pad = 32

    out = Matrix{Float32}(undef, dims, N)
    bitarray = Array{Int,3}(undef, pad, N, dims)
    sobols!(rng, out, bitarray)

    return out
end

sobols(dims; m=9) = sobols(Random.Xoshiro(31), dims; m=m)

"""
Sample a Sobol' sequence into a preallocated matrix and scramble it.
"""
function sobols!(
    rng::Random.AbstractRNG,
    out::AbstractMatrix{T},
    bitarray::AbstractArray{<:Integer,3},
) where {T<:Real}
    dims, N = size(out)

    try
        m = Int(log2(N))
    catch err
        if err isa InexactError
            throw(
                ArgumentError("Number of matrix columns is not a power of 2"),
            )
        else
            rethrow()
        end
    end

    if N == 0
        return out
    end

    _sobols!(out)
    _scramble!(rng, out, bitarray)

    return out
end

"""
Sample a Sobol' sequence into a preallocated matrix.
"""
function _sobols!(out::AbstractMatrix{T}) where {T<:AbstractFloat}
    dims, N = size(out)
    seq = Sobol.SobolSeq(dims)
    # `Sobol.jl` does not include `(0)^d` but this invalidates the digital net
    # property. The common practice to work around that is to draw `2N` and
    # throw away the first `N` Sobol' points. While this is OK for low `N`, it
    # is overly wasteful for larger `N`. This is addressed as well by this (as
    # of 2025-03-21 open) PR: https://github.com/JuliaMath/Sobol.jl/pull/35 .
    #
    # See https://github.com/JuliaMath/Sobol.jl/issues/31#issue-1328916662
    # (QuasiMonteCarlo.jl is based on Sobol.jl) and
    # https://github.com/JuliaMath/Sobol.jl/blob/685cec3fde77dce494c208f2de36c89f438254f6/src/Sobol.jl#L77
    #
    # My solution here is to simply do what the PR does and include `(0)^d`.
    # Since the sequence is scrambled anyways, this should not introduce any
    # weirdness due to the zeroes (to cite from previous links, “But as Owen's
    # paper shows (Figure 1), an even better solution is to scramble the
    # points.”).
    for (idx, x) in enumerate(eachcol(out))
        if idx == 1
            x .= zeros(length(x))
        else
            Sobol.next!(seq, x)
        end
    end
    return out
end

"""
Scramble a Sobol' sequence using preallocated matrices.
"""
function _scramble!(
    rng::Random.AbstractRNG,
    out::AbstractMatrix{T},
    bitarray::AbstractArray{<:Integer,3},
) where {T<:Number}
    dims, N = size(out)

    # This is the default for the Matoušek scramble.
    pad = 32
    # This is fixed for the Matoušek scramble of the Sobol' sequence.
    b = 2

    @assert size(bitarray) == (pad, N, dims)

    # Note that `bitarray` has swapped `N` and `dims` due to how it is later
    # indexed into (we want to iterate over the first dimensions first). This
    # means that we have to transpose `out` which is not as efficient as it may
    # be able to be.
    out_ = out'

    QuasiMonteCarlo.unif2bits!(bitarray, out_, b)
    randomize_bits!(bitarray, MatousekScramble(; base=b, pad=pad, rng=rng))
    for i in CartesianIndices(out_)
        out_[i] = QuasiMonteCarlo.bits2unif(T, @view(bitarray[:, i]), b)
    end

    return out
end

"""
Truly in-place version of `randomize_bits!` (i.e. overwrites the original array).

Copied and adjusted from
[QuasiMonteCarlo.jl](https://github.com/SciML/QuasiMonteCarlo.jl).
"""
function randomize_bits!(
    bits::AbstractArray{T,3},
    R::MatousekScramble,
) where {T<:Integer}
    # https://statweb.stanford.edu/~owen/mc/ Chapter 17.6 around equation (17.15).
    #
    pad, n, d = size(bits)
    b = R.base
    rng = R.rng
    m = QuasiMonteCarlo.logi(b, n)
    @assert m ≥ 1 "We need m ≥ 1" # m=0 causes awkward corner case below.  Caller handles that case specially.

    for s in 1:d
        # Permutations matrix and shift to apply to bits 1:m
        matousek_M, matousek_C = QuasiMonteCarlo.getmatousek(rng, m, b)

        # xₖ = (∑ₗ Mₖₗ aₗ + Cₖ) mod b where xₖ is the k element in base b
        # matousek_M (m×m) * origin_bits (m×n) .+ matousek_C (m×1)
        @views bits[1:m, :, s] .=
            (matousek_M * bits[1:m, :, s] .+ matousek_C) .% b
    end

    # Paste in random entries for bits after m'th one
    if pad > m
        # random_bits[(m + 1):pad, :, :] = rand(rng, 0:(b - 1), n * d * (pad - m))
        rand!(rng, @view(bits[(m + 1):pad, :, :]), 0:(b - 1))
    end
end
