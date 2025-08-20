using LandscapeSimple
using Random

# Note that we use a `NamedTuple`s instead of `Dict`s to keep the order.
hps_vary = (;
    a=mkscale_discrete([1.0, 2.0, 4.0]),
    b=mkscale_minmax(1.0, 5.0),
    c=mkscale_geo(0.0, 3.0),
    d=ceil ∘ mkscale_geo(0.0, 4.0),
    e=mkscale_mix([0.2 => mkscale_const(0.0), 0.8 => mkscale_geo(0.0, 5.0)]),
)

# Generate `2^m` configurations deterministically.
m = 5
configs = configurations(Random.Xoshiro(31), hps_vary; m=m)

try
    using CairoMakie

    table = Matrix{Real}(undef, length(configs), length(hps_vary))
    for (idx_name, name) in enumerate(propertynames(hps_vary))
        table[:, idx_name] .= getproperty.(configs, Ref(name))
    end

    labels = propertynames(hps_vary)
    idxs_label = repeat(1:length(labels); inner=length(configs))
    fig = Figure()
    ax = Axis(fig[1, 1]; xticks=(1:length(labels), string.(collect(labels))))
    scatter!(
        ax,
        idxs_label .+ randn(length(idxs_label)) .* 0.01,
        table[:];
        marker='+',
    )
    display(fig)

catch e
    if e isa ArgumentError
        @info "CairoMakie not available, skipping plots"
    else
        rethrow()
    end
end
