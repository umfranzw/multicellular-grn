using IndividualMod
using RunMod

run = RunMod.get_first_run()
pop = Array{Individual, 1}()
for i in 1:2
    indiv = IndividualMod.rand_init(run, UInt64(i))
    push!(pop, indiv)
end

prev_pop = deepcopy(pop)
pop[1].genes[1].bind_sites[1].tag += 1

println(prev_pop[1].genes[1].bind_sites[1].tag)
println(pop[1].genes[1].bind_sites[1].tag)
