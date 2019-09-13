import Random

rand_bits(run::RunInfo.Run, n::Int64) = Random.bitrand(run.rng, n)
rand_int(run::RunInfo.Run, low::Int64, high::Int64) = Random.rand(run.rng, low:high)
rand_float(run::RunInfo.Run) = Random.rand(run.rng)
