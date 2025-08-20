using LandscapeSimple
using Random
using StatsBase
using Test

@testset "Configuration Generator Tests" begin
    @testset "Basic Configuration Generation" begin
        # Setup hyperparameter scales.
        hps_vary = (;
            a=mkscale_discrete([3, 5, 8, 12, 17, 23]),
            c=mkscale_geo(0.0, 2.0),
            d=Int ∘ ceil ∘ mkscale_geo(0.0f0, 10.0f0),
            b=mkscale_minmax(1.0e-3, 5.0e-2),
        )

        # Generate configurations.
        m = 7
        configs = configurations(Random.Xoshiro(31), hps_vary; m=m)

        @test length(configs) == 2^m
        @test length(configs) == 128

        # Test that all configurations have the expected fields.
        @test all(hasfield(typeof(c), :a) for c in configs)
        @test all(hasfield(typeof(c), :b) for c in configs)
        @test all(hasfield(typeof(c), :c) for c in configs)
        @test all(hasfield(typeof(c), :d) for c in configs)
    end

    @testset "Deterministic Behavior" begin
        hps_vary = (;
            a=mkscale_discrete([3, 5, 8, 12, 17, 23]),
            c=mkscale_geo(0.0, 2.0),
            d=Int ∘ ceil ∘ mkscale_geo(0.0f0, 10.0f0),
            b=mkscale_minmax(1.0e-3, 5.0e-2),
        )

        # Generate configurations twice with same seed.
        configs1 = configurations(Random.Xoshiro(31), hps_vary; m=5)
        configs2 = configurations(Random.Xoshiro(31), hps_vary; m=5)

        # Should produce identical results.
        @test length(configs1) == length(configs2)
        for i in 1:length(configs1)
            @test configs1[i].a == configs2[i].a
            @test configs1[i].b == configs2[i].b
            @test configs1[i].c == configs2[i].c
            @test configs1[i].d == configs2[i].d
        end

        # Different seed should produce different results.
        configs3 = configurations(Random.Xoshiro(42), hps_vary; m=5)
        @test any(configs1[i] != configs3[i] for i in 1:length(configs1))
    end

    @testset "Value Range Validation" begin
        valid_a_values = [3, 5, 8, 12, 17, 23]
        lc, uc = 0.0, 2.0
        lb, ub = 1.0e-3, 5.0e-2
        ld, ud = 1.0, 10.0
        hps_vary = (;
            a=mkscale_discrete(valid_a_values),
            c=mkscale_geo(lc, uc),
            d=Int ∘ ceil ∘ mkscale_geo(ld, ud),
            b=mkscale_minmax(lb, ub),
        )

        configs = configurations(Random.Xoshiro(31), hps_vary; m=8)

        # Test discrete values for 'a'.
        @test all(c.a in valid_a_values for c in configs)

        # Test range for 'b' (minmax scale).
        @test all(lb <= c.b <= ub for c in configs)

        # Test range for 'c' (geometric scale).
        @test all(lc <= c.c <= uc for c in configs)

        # Test 'd' is integer and in valid range.
        @test all(isa(c.d, Int) for c in configs)
        @test all(1 <= c.d <= 10 for c in configs)
    end

    @testset "Edge Cases" begin
        # TODO Add geo edge cases (0, negative values, …)

        # Test with single value discrete scale.
        hps_single = (;
            a=mkscale_discrete([42]),
            b=mkscale_minmax(1.0, 1.0),  # Single value range
        )

        configs_single = configurations(Random.Xoshiro(31), hps_single; m=3)
        @test all(c.a == 42 for c in configs_single)
        @test all(c.b == 1.0 for c in configs_single)
    end

    @testset "Distribution Properties" begin
        hps_vary = (;
            a=mkscale_discrete([1, 2, 3, 4, 5, 6]),
            b=mkscale_minmax(0.0, 1.0),
        )

        # Generate many configurations.
        configs = configurations(Random.Xoshiro(31), hps_vary; m=10)

        # Check that discrete values are reasonably distributed.
        a_counts = Dict(i => count(c.a == i for c in configs) for i in 1:6)

        # Each value should appear at least once in 1024 samples.
        @test all(count > 0 for count in values(a_counts))

        # Check continuous values have reasonable spread.
        b_values = [c.b for c in configs]
        @test minimum(b_values) < 0.1  # Some values near bottom of range
        @test maximum(b_values) > 0.9  # Some values near top of range
        @test 0.4 < mean(b_values) < 0.6  # Mean roughly centered
    end

    @testset "Type Stability" begin
        hps_vary = (;
            a=mkscale_discrete([3, 5, 8, 12, 17, 23]),
            c=mkscale_geo(0.0, 2.0),
            d=Int ∘ ceil ∘ mkscale_geo(0.0f0, 10.0f0),
            b=mkscale_minmax(1.0e-3, 5.0e-2),
        )

        configs = configurations(Random.Xoshiro(31), hps_vary; m=5)

        # Check that all configs have consistent types.
        first_config = configs[1]
        for config in configs
            @test typeof(config.a) == typeof(first_config.a)
            @test typeof(config.b) == typeof(first_config.b)
            @test typeof(config.c) == typeof(first_config.c)
            @test typeof(config.d) == typeof(first_config.d)
        end

        # Verify specific expected types.
        @test isa(configs[1].a, Integer)
        @test isa(configs[1].b, AbstractFloat)
        @test isa(configs[1].c, AbstractFloat)
        @test isa(configs[1].d, Int)
    end

    @testset "Output Type Preservation" begin

        # Test with Int32 discrete values.
        hps_int32 = (; a=mkscale_discrete(Int32[10, 20, 30, 40]),)
        configs_int32 = configurations(Random.Xoshiro(42), hps_int32; m=4)
        @test all(isa(c.a, Int32) for c in configs_int32)
        @test eltype([c.a for c in configs_int32]) == Int32

        # Test with Int64 discrete values.
        hps_int64 = (; a=mkscale_discrete(Int64[100, 200, 300]),)
        configs_int64 = configurations(Random.Xoshiro(42), hps_int64; m=4)
        @test all(isa(c.a, Int64) for c in configs_int64)
        @test eltype([c.a for c in configs_int64]) == Int64

        # Test with Float32 values.
        hps_float32 = (;
            b=mkscale_minmax(Float32(0.1), Float32(0.9)),
            c=mkscale_geo(Float32(1.0), Float32(10.0)),
        )
        configs_float32 = configurations(Random.Xoshiro(42), hps_float32; m=4)
        @test all(isa(c.b, Float32) for c in configs_float32)
        @test all(isa(c.c, Float32) for c in configs_float32)
        @test eltype([c.b for c in configs_float32]) == Float32
        @test eltype([c.c for c in configs_float32]) == Float32

        # Test with Float64 values.
        hps_float64 = (;
            b=mkscale_minmax(0.1, 0.9),  # Float64 by default
            c=mkscale_geo(1.0, 10.0),    # Float64 by default
        )
        configs_float64 = configurations(Random.Xoshiro(42), hps_float64; m=4)
        @test all(isa(c.b, Float64) for c in configs_float64)
        @test all(isa(c.c, Float64) for c in configs_float64)
        @test eltype([c.b for c in configs_float64]) == Float64
        @test eltype([c.c for c in configs_float64]) == Float64

        # Test with mixed numeric types.
        hps_mixed = (;
            a=mkscale_discrete(UInt8[1, 2, 3, 4]),
            b=mkscale_discrete(Int16[10, 20, 30]),
            c=mkscale_minmax(Float16(0.0), Float16(1.0)),
            d=mkscale_geo(BigFloat(1.0), BigFloat(100.0)),
        )
        configs_mixed = configurations(Random.Xoshiro(42), hps_mixed; m=3)

        @test all(isa(c.a, UInt8) for c in configs_mixed)
        @test all(isa(c.b, Int16) for c in configs_mixed)
        @test all(isa(c.c, Float16) for c in configs_mixed)
        @test all(isa(c.d, BigFloat) for c in configs_mixed)

        # Test with custom types.
        hps_custom = (;
            symbols=mkscale_discrete([:low, :medium, :high]),
            strings=mkscale_discrete(["option1", "option2", "option3"]),
        )
        configs_custom = configurations(Random.Xoshiro(42), hps_custom; m=3)

        @test all(isa(c.symbols, Symbol) for c in configs_custom)
        @test all(isa(c.strings, String) for c in configs_custom)

        # Test composed functions preserve types
        hps_composed = (;
            # Float32 input should result in Int32 output after composition
            a=Int32 ∘ round ∘ mkscale_minmax(Float32(0.0), Float32(10.0)),
            # Float64 input should result in Int64 output after composition
            b=Int64 ∘ floor ∘ mkscale_geo(1.0, 100.0),
            # Ensure ceil preserves input type
            c=ceil ∘ mkscale_minmax(Float32(0.0), Float32(5.0)),
        )
        configs_composed =
            configurations(Random.Xoshiro(42), hps_composed; m=4)

        @test all(isa(c.a, Int32) for c in configs_composed)
        @test all(isa(c.b, Int64) for c in configs_composed)
        @test all(isa(c.c, Float32) for c in configs_composed)
    end

    @testset "Large Scale Generation" begin
        hps_vary =
            (; a=mkscale_discrete(collect(1:100)), b=mkscale_minmax(0.0, 1.0))

        # Test with larger m value.
        m = 12  # 4096 configurations
        configs = configurations(Random.Xoshiro(31), hps_vary; m=m)

        @test length(configs) == 2^m
        @test length(configs) == 4096

        # Verify memory efficiency - all configs should be unique.
        unique_configs = unique(configs)
        @test length(unique_configs) == length(configs)
    end
end

# --- Minimal mock "scales" for testing ---------------------------------------

# Callable scale: returns a tagged value so we can see which scale was used.
callable_scale(tag) = u -> (tag, u)  # returns (which_scale, local_u)

# Object with a `.transform` method:
struct ObjScale
    tag::Symbol
    transform
end
ObjScale(tag) = ObjScale(tag, x -> (tag, x))
(transform::ObjScale)(u) = error("should not call as a function")

# Utility to extract tag from mix output (first element of tuple)
get_tag(x) = x[1]
get_u(x) = x[2]

# --- Tests --------------------------------------------------------------------

@testset "mkscale_mix" begin
    # 1) Basic selection by proportions, with normalisation
    @testset "bin selection & normalisation" begin
        scales = [
            0.2 => callable_scale(:A),
            0.5 => ObjScale(:B),
            0.3 => callable_scale(:C),
        ]
        mix = mkscale_mix(scales; normalize=true)

        # Representative points inside each bin (using p cumulative: 0.2, 0.7, 1.0)
        @test get_tag(mix(0.00)) === :A
        @test get_tag(mix(0.199999999)) === :A  # just below first edge
        @test get_tag(mix(0.200000001)) === :B
        @test get_tag(mix(0.699999999)) === :B
        @test get_tag(mix(0.700000001)) === :C
        @test get_tag(mix(0.999999999)) === :C

        # Right edge inclusivity: last bin includes u == 1.0
        @test get_tag(mix(1.0)) === :C
    end

    # 2) Within-bin local rescaling is correct (v = (u - l)/(u - l))
    @testset "within-bin rescaling" begin
        scales = [0.5 => callable_scale(:L), 0.5 => callable_scale(:R)]
        mix = mkscale_mix(scales)  # edges: [0.5, 1.0]

        # Midpoints of each bin should map to v ≈ 0.5
        tagL, vL = mix(0.25)
        @test tagL === :L
        @test isapprox(vL, 0.5; atol=1e-12)

        tagR, vR = mix(0.75)
        @test tagR === :R
        @test isapprox(vR, 0.5; atol=1e-12)

        # Left edge maps to v = 0, right bin left edge maps to v = 0
        _, v0 = mix(0.0)
        @test isapprox(v0, 0.0; atol=1e-12)
        _, vR0 = mix(0.5)
        @test isapprox(vR0, 0.0; atol=1e-12)

        # Values extremely close to 1 map to last bin with v < 1
        _, vTop = mix(1.0 - eps(Float64))
        @test vTop < 1.0
    end

    # 3) Callable vs .transform dispatch both work
    @testset "callable and .transform both supported" begin
        scales = [0.3 => callable_scale(:F), 0.7 => ObjScale(:T)]
        mix = mkscale_mix(scales)

        @test get_tag(mix(0.1)) === :F    # first component (callable)
        @test get_tag(mix(0.9)) === :T    # second component (.transform)
    end

    # 4) Normalise=false requires sum≈1, otherwise error
    @testset "normalise flag / sum≈1 enforcement" begin
        scales_ok = [0.25 => callable_scale(:X), 0.75 => callable_scale(:Y)]
        @testset "sum exactly 1" begin
            mix_ok = mkscale_mix(scales_ok; normalize=false)
            @test get_tag(mix_ok(0.2)) === :X
            @test get_tag(mix_ok(0.8)) === :Y
        end

        scales_bad = [0.25 => callable_scale(:X), 0.80 => callable_scale(:Y)] # sum=1.05
        @test_throws AssertionError mkscale_mix(scales_bad; normalize=false)
    end

    # 5) Proportion validation (non-positive and non-finite)
    @testset "invalid proportions" begin
        @test_throws AssertionError mkscale_mix([0.0 => callable_scale(:A)])
        @test_throws AssertionError mkscale_mix([
            -0.1 => callable_scale(:A),
            1.1 => callable_scale(:B),
        ])
        @test_throws AssertionError mkscale_mix([
            NaN => callable_scale(:A),
            1.0 => callable_scale(:B),
        ])
    end

    # 6) Approximate frequency check over a uniform grid (not Sobol’ specific)
    @testset "empirical proportions over grid" begin
        scales = [
            0.2 => callable_scale(:A),
            0.5 => callable_scale(:B),
            0.3 => callable_scale(:C),
        ]
        mix = mkscale_mix(scales)

        N = 10_000
        counts = Dict(:A => 0, :B => 0, :C => 0)
        for k in 0:(N - 1)
            u = (k + 0.5) / N  # midpoints
            counts[get_tag(mix(u))] += 1
        end
        # Expect close to target proportions within a small tolerance
        @test isapprox(counts[:A] / N, 0.2; atol=5e-3)
        @test isapprox(counts[:B] / N, 0.5; atol=5e-3)
        @test isapprox(counts[:C] / N, 0.3; atol=5e-3)
    end

    # 7) Boundary and clamping behaviour
    @testset "boundary handling" begin
        scales = [0.9 => callable_scale(:M), 0.1 => callable_scale(:N)]
        mix = mkscale_mix(scales)

        @test get_tag(mix(0.0)) === :M
        @test get_tag(mix(nextfloat(0.9))) === :N   # just over the first edge
        @test get_tag(mix(1.0)) === :N              # exact 1.0 -> last bin
    end
end
