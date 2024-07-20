using LandscapeSimple
using Random

# Note that we use a `NamedTuple`s instead of `Dict`s to keep the order.
hps_vary = (;
    a=mkscale_discrete([3, 5, 8, 12, 17, 23]),
    b=mkscale_minmax(1.0f-3, 5.0f-2),
)

# Generate `2^m` configurations deterministically. For now, we never run
# more than 1024 configs.
m = 5
configs = configurations(Random.Xoshiro(31), hps_vary; m=m)

# If you have UnicodePlots available:
try
    using UnicodePlots
    local a = map(first, configs)
    local b = map(last, configs)
    UnicodePlots.scatterplot(a, b)
catch
    println("Skipping UnicodePlots stuff")
end
