module RunMod

import TOML
import Random

export Run,
    get_config_channel, get_first_run

const CONFIG_PATH = "/home/umfranzw/multicellular-grn/runs.toml"

struct Run
    pop_size::Int64
    ea_steps::Int64
    mut_prob::Float64
    num_genes::Int64
    num_bind_sites::Int64
    fitness_term_threshold::Float64

    reg_steps::Int64
    min_protein_conc::Float64
    max_protein_conc::Float64
    output_rate_step::Float64
    genes_per_cell::Int64

    decay_rate::Float64
    num_initial_tfs::Int64
    max_proteins::Int64
    max_mut_delta::Float64
    max_mut_bits::Int64
    binding_seq_play::Int64

    growth_start::Int64
    growth_end::Int64
    growth_threshold::Float64
    diff_threshold::Float64

    code_start::Int64
    code_end::Int64

    initial_cell_energy::Float64

    fix_rng_seed::Bool
    rng_seed::Int64

    log_level::String
    log_buff_size::Int64
    log_dir::String

    rng::Random.MersenneTwister

    function Run(run)
        new(
            run["pop_size"],
            run["ea_steps"],
            run["mut_prob"],
            run["num_genes"],
            run["num_bind_sites"],
            run["fitness_term_threshold"],
            
            run["reg_steps"],
            run["min_protein_conc"],
            run["max_protein_conc"],
            run["output_rate_step"],
            run["genes_per_cell"],
            
            run["decay_rate"],
            run["num_initial_tfs"],
            run["max_proteins"],
            run["max_mut_delta"],
            run["max_mut_bits"],
            
            run["binding_seq_play"],
            
            run["growth_start"],
            run["growth_end"],
            run["growth_threshold"],
            run["diff_threshold"],
            
            run["code_start"],
            run["code_end"],
            
            run["initial_cell_energy"],
            
            run["fix_rng_seed"],
            run["rng_seed"],
            
            run["log_level"],
            run["log_buff_size"],
            run["log_dir"],
            
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
