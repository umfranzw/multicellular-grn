module RandUtilsMod

using RunMod
import Random

export rand_bits, rand_int, rand_float, rand_floats, rand_enum_val

rand_bits(config::Config, n::Int64) = Random.bitrand(config.rng, n)
rand_int(config::Config, low::Int64, high::Int64) = Random.rand(config.rng, low:high) #range is [low, high]
rand_float(config::Config) = Random.rand(config.rng) #range is [0, 1)
rand_floats(config::Config, n::Int64) = Random.rand(config.rng, n)
rand_floats(config::Config, start::Float64, stop::Float64, n::Int64) = map(r -> start + r * (stop - start), Random.rand(config.rng, n))
rand_enum_val(config::Config, enum::DataType) = Random.rand(config.rng, instances(enum))

end
