module RunMod

import TOML
import Random

export Run,
    get_config_channel, get_first_run

#const CONFIG_PATH = "/home/umfranzw/multicellular-grn/runs.toml"
const CONFIG_PATH = "/home/wayne/Documents/school/thesis/multicellular-grn/runs.toml"

struct Run
    pop_size::Int64
    ea_steps::Int64
    mut_prob::Float64
    num_genes::Int64
    num_bind_sites::Int64
    fitness_term_threshold::Float64
    reg_bind_threshold::Float64
    growth_bind_threshold::Float64
    bind_bind_threshold::Float64
    prod_bind_threshold::Float64
    prod_rate_incr::Float64

    reg_steps::Int64
    min_protein_threshold::Float64

    decay_rate::Float64
    num_initial_proteins::Int64
    max_proteins::Int64
    max_mut_float_change::Float64
    max_mut_bits_change::Int64
    binding_seq_play::Int64

    diff_alpha::Float64
    diff_h::Float64
    diff_dt::Float64

    dev_start_iter::Int64
    cell_div_conc_threshold::Float64

    initial_cell_energy::Float64

    fix_rng_seed::Bool
    rng_seed::Int64

    rng::Random.MersenneTwister

    function Run(run)
        new(
            run["pop_size"],
            run["ea_steps"],
            run["mut_prob"],
            run["num_genes"],
            run["num_bind_sites"],
            run["fitness_term_threshold"],
            run["reg_bind_threshold"],
            run["growth_bind_threshold"],
            run["bind_bind_threshold"],
            run["prod_bind_threshold"],
            run["prod_rate_incr"],
            
            run["reg_steps"],
            run["min_protein_threshold"],
            
            run["decay_rate"],
            run["num_initial_proteins"],
            run["max_proteins"],
            run["max_mut_float_change"],
            run["max_mut_bits_change"],
            run["binding_seq_play"],

            run["diff_alpha"],
            run["diff_h"],
            run["diff_dt"],
            
            run["dev_start_iter"],
            run["cell_div_conc_threshold"],
            
            run["initial_cell_energy"],
            
            run["fix_rng_seed"],
            run["rng_seed"],
            
            Random.MersenneTwister()
        )
    end
end

function next_run(c::Channel)
    config = TOML.parsefile(CONFIG_PATH)
    for i in 1:length(config["runs"])
        put!(c, Run(config["runs"][i]))
    end
end


function get_config_channel()
    Channel(next_run)
end

#useful for testing in the REPL
function get_first_run()
    channel = get_config_channel()
    take!(channel)
end

end
