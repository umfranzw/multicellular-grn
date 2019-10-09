import RunMod
#using IndividualMod
using ProteinPropsMod

#run = RunMod.get_first_run()
#indiv = IndividualMod.rand_init(run)
p = ProteinProps(ProteinPropsMod.Reg, ProteinPropsMod.Inter, ProteinPropsMod.Activate, ProteinPropsMod.D)
println(p)
