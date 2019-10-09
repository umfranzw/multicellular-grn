import RunMod
using IndividualMod

run = RunMod.get_first_run()
indiv = IndividualMod.rand_init(run)
println(indiv)
