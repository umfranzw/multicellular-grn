module RunMod

import TOML
import Random

export Run, Config,
    get_run_channel, get_first_run

if gethostname() == "ibis"
    const CONFIG_PATH = "/home/umfranzw/multicellular-grn/runs.toml"
    #no trailing "/"
    const DATA_PATH = "/home/umfranzw/multicellular-grn/data"
else
    const CONFIG_PATH = "/home/wayne/Documents/school/thesis/multicellular-grn/runs.toml"
    #no trailing "/"
    const DATA_PATH = "/home/wayne/Documents/school/thesis/multicellular-grn/data"
end

struct Run
    pop_size::Int64
    ea_steps::Int64
    mut_prob::Float64
    num_initial_genes::Int64
    fitness_term_threshold::Float64
    reg_bind_threshold::Float64
    protein_app_threshold::Float64
    prod_rate_incr::Float64

    reg_steps::Int64
    protein_deletion_threshold::Float64

    decay_rate::Float64
    num_initial_proteins::Int64
    max_proteins::Int64
    max_conc_mut::Float64
    max_mut_bits_change::Int64

    diff_alpha::Float64
    diff_h::Float64
    diff_dt::Float64

    accuracy_weight::Float64
    evolvability_weight::Float64

    cell_energy_threshold::Float64

    initial_cell_energy::Float64

    fix_rng_seed::Bool
    rng_seed::UInt64

    log_data::Bool
    step_range::StepRange{Int64, Int64}
    data_output_file::String

    function Run(run)
        new(
            run["pop_size"],
            run["ea_steps"],
            run["mut_prob"],
            run["num_initial_genes"],
            run["fitness_term_threshold"],
            run["reg_bind_threshold"],
            run["protein_app_threshold"],
            run["prod_rate_incr"],
            
            run["reg_steps"],
            run["protein_deletion_threshold"],
            
            run["decay_rate"],
            run["num_initial_proteins"],
            run["max_proteins"],
            run["max_conc_mut"],
            run["max_mut_bits_change"],

            run["diff_alpha"],
            run["diff_h"],
            run["diff_dt"],

            run["accuracy_weight"],
            run["evolvability_weight"],
            
            run["cell_energy_threshold"],
            
            run["initial_cell_energy"],
            
            run["fix_rng_seed"],
            run["rng_seed"],

            run["log_data"],
            parse_step_range(run["step_range"]),
            run["data_output_file"]
        )
    end
end

mutable struct Config
    run::Run
    rng::Random.MersenneTwister
end

function parse_step_range(str::String)
    #note: could just eval the string, but this allows us to give a more descriptive error message when something goes wrong
    try
        vals = map(val -> parse(Int64, val), split(str, ":"))
        return StepRange(vals...)
    catch
        @error("Error parsing step_range param from run file.")
        exit(1)
    end
end

function next_run(c::Channel)
    config = TOML.parsefile(CONFIG_PATH)
    for i in 1:length(config["runs"])
        put!(c, Run(config["runs"][i]))
    end
end


function get_run_channel()
    Channel(next_run)
end

#useful for testing in the REPL
function get_first_run()
    channel = get_run_channel()
    take!(channel)
end

function get_test_config()
    run = get_first_run()

    Config(run, Random.MersenneTwister())
end

end
