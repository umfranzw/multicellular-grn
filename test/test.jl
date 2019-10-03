import RunMod
using ProteinMod

run = RunMod.get_first_run()
p = Protein(run, ProteinProps(ProteinMod.Reg, ProteinMod.Intra, ProteinMod.RateUp, ProteinMod.A), true)
