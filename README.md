# LandscapeSimple.jl


A tiny simplistic hackable library that generates lists of configurations of a
mix of continuous and discrete parameters using Quasi–Monte Carlo methods to
cover the space faster than uniform random sampling. For
now, based on the Sobol' sequence. Based on
[QuasiMonteCarlo.jl](https://github.com/SciML/QuasiMonteCarlo.jl)'s Sobol'
sequence which is in turn based on
[Sobol.jl](https://github.com/JuliaMath/Sobol.jl).  Uses the following defaults,
which seem to be sane (please open an issue if you disagree):

- Skip the first $2^n-1$ in the sequence due to `Sobol.jl` not starting at
  $(0)^d$ (see also [this
  issue](https://github.com/JuliaMath/Sobol.jl/issues/31#issue-1328916662)).
- Perform sequence scrambling via `QuasiMonteCarlo.jl`.

Example (copy-pasted from `scripts/example.jl` for convenience):

```Julia
using LandscapeSimple
using Random

# Note that we use a `NamedTuple`s instead of `Dict`s to keep the order.
hps_vary = (;
    a=mkscale_discrete([3, 5, 8, 12, 17, 23]),
    c=mkscale_geo(0.0f0, 2.0f0),
    d=Int ∘ ceil ∘ mkscale_geo(0.0f0, 2.0f0),
    b=mkscale_minmax(1.0f-3, 5.0f-2),
)

# Generate `2^m` configurations deterministically.
m = 7
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
        println(configs)
    else
        rethrow()
    end
end
```

yields

```
Uniformly sampled tends to clump:
        ┌────────────────────────────────────────┐ 
   0.05 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠠⠀⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
        │⠀⠀⠀⠀⡀⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠈⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠐⠀⠀⠀⢀⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⢀⠀⠀│ 
        │⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠀⠀│ 
        │⠀⠀⠀⠀⠄⠀⠀⠀⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠄⠀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠠⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠰⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
      0 │⠀⠀⠀⠀⠄⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
        └────────────────────────────────────────┘ 
        ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀30⠀ 
Sobol' sequence more evenly covers space:
        ┌────────────────────────────────────────┐ 
   0.05 │⠀⠀⠀⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⡀⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
        │⠀⠀⠀⠀⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠨⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⠀⠀│ 
        │⠀⠀⠀⠀⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
        │⠀⠀⠀⠀⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
        │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
        │⠀⠀⠀⠀⡀⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⠀⠀│ 
      0 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
        └────────────────────────────────────────┘ 
        ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀30⠀ 
```


## Warning


The Sobol' sequence has nice properties but only if you adhere to a few rules.
Most importantly, you have to use *all of the first $2^n$ elements* for the
nicest properties. This seems to be especially important for discrete spaces.
Read about the Sobol' sequence before relying on this code.


## Open points/TODOs


- [ ] Allow to provide `Float64` to `mkscale_minmax` to get `Float64` out

  Right now the underlying Sobol' sequence is hardcoded to `Float32` since what
  it was originally intended to be used for was `Float32`. This means that the
  numbers resulting from `mkscale_minmax` (as well as `mkscale_geo` but for that
  function this issue is currently not made explicit) only can have a `Float32`
  resolution which may not always be what you want.
