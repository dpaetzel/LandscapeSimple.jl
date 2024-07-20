using LandscapeSimple
using Random

# Note that we use a `NamedTuple`s instead of `Dict`s to keep the order.
hps_vary = (;
    a=mkscale_discrete([3, 5, 8, 12, 17, 23]),
    b=mkscale_minmax(1.0f-3, 5.0f-2),
)

# Generate `2^m` configurations deterministically.
m = 5
configs = configurations(Random.Xoshiro(31), hps_vary; m=m)

# If you have UnicodePlots available:
try
    using StatsBase
    using UnicodePlots

    println("Uniformly sampled tends to clump:")
    a_ = sample([3, 5, 8, 12, 17, 23], 2^m)
    b_ = rand(2^m) .* (5.0f-2 - 1.0f-3) .+ 1.0f-3
    display(UnicodePlots.scatterplot(Float64.(a_), b_))

    println("Sobol' sequence more evenly covers space:")
    a = map(first, configs)
    b = map(last, configs)
    display(UnicodePlots.scatterplot(a, b))
catch e
    if e isa ArgumentError
        println("Skipping StatsBase/UnicodePlots stuff")
    else
        rethrow()
    end
end
