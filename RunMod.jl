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
    max_genes::Int64
    tourn_size::Int64
    fitness_term_threshold::Float64
    gene_score_threshold::Float64
    
    reg_steps::Int64
    protein_deletion_threshold::Float64
    max_prod_rate::Float64

    decay_rate::Float64
    num_initial_proteins::Int64
    initial_cell_sensor_conc::Float64
    max_proteins_per_cell::Int64
    
    diff_alpha::Float64
    diff_h::Float64
    diff_dt::Float64

    initial_acc_weight::Float64
    initial_ev_weight::Float64
    weight_shift::Float64
    max_acc_weight::Float64
    min_ev_weight::Float64

    bind_sites_per_gene::Int64

    division_age_limit::Int64
    cell_division_threshold::Float64
    max_children::Int64
    max_tree_size::Int64
    max_sensor_amount::Float64
    sensor_reinforcement_threshold::Float64
    sym_prob_threshold::Float64
    max_protein_arg::UInt8
    
    fix_rng_seed::Bool
    rng_seed::UInt64

    log_data::Bool
    step_range::StepRange{Int64, Int64}
    data_output_file::String

    function Run(run::Dict{AbstractString, Any})
        args = Array{Any, 1}()
        for name_sym in fieldnames(Run)
            name_str = string(name_sym)
            if name_str == "step_range"
                push!(args, parse_step_range(run[name_str]))
            else
                push!(args, run[name_str])
            end
        end
        
        new(args...)
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
