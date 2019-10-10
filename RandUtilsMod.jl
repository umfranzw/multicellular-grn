module RandUtilsMod

using RunMod
import Random

export rand_bits, rand_int, rand_float, rand_floats, rand_enum_val

rand_bits(run::Run, n::Int64) = Random.bitrand(run.rng, n)
rand_int(run::Run, low::Int64, high::Int64) = Random.rand(run.rng, low:high) #range is [low, high]
rand_float(run::Run) = Random.rand(run.rng) #range is [0, 1)
rand_floats(run::Run, n::Int64) = Random.rand(run.rng, n)
rand_enum_val(run::Run, enum::DataType) = Random.rand(run.rng, instances(enum))

end
