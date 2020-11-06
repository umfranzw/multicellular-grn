module RunMod

import TOML
import Random
import CompressionMod

export Run, Config,
    get_run_channel, get_first_run

if gethostname() == "ibis"
    const CONFIG_PATH = "/home/umfranzw/multicellular-grn/runs.toml"
    #no trailing "/"
    const DATA_PATH = "/var/www/data"
else
    const CONFIG_PATH = "/home/wayne/Documents/school/thesis/multicellular-grn/runs.toml"
    #no trailing "/"
    const DATA_PATH = "/home/wayne/Documents/school/thesis/multicellular-grn/data"
end

@enum LogLevel::UInt8 LogNone=0 LogFitnesses LogIndivs

struct Run
    pop_size::Int64
    ea_steps::Int64
    point_mut_prob::Float64
    dup_mut_prob::Float64
    sys_level_mut_prob::Float64
    cross_prop::Float64
    growth_threshold::Float64
    num_initial_genes::Int64
    max_genes::Int64
    tourn_size::Int64
    fitness_term_threshold::Float64
    gene_dup_count_threshold::Int64
    gene_sys_level_threshold::Int64
    gene_dup_gen_limit::Int64
    gene_sys_level_gen_limit::Int64

    tag_limit::Int64
    
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

    bind_threshold::Float64
    bind_consum_rate::Float64
        
    division_age_limit::Int64
    cell_division_threshold::Float64
    max_children::Int64
    max_tree_size::Int64
    max_sensor_amount::Float64
    sensor_reinforcement_threshold::Float64
    sym_prob_threshold::Float64
    max_protein_arg::UInt8
    
    fix_rng_seed::Bool
    user_rng_seed::UInt64

    log_level::LogLevel
    step_range::StepRange{Int64, Int64}
    data_output_file::String
    multithreaded::Bool
    compression_alg::CompressionMod.CompressionAlg

    #note: different TOML versions use AbstractString or String
    function Run(run::Union{Dict{AbstractString, Any}, Dict{String, Any}})
        compression_alg_dict = get_enum_dict(CompressionMod.CompressionAlg)
        
        args = Array{Any, 1}()
        for name_sym in fieldnames(Run)
            name_str = string(name_sym)
            if name_str == "step_range"
                push!(args, parse_step_range(run[name_str]))
            elseif name_str == "log_level"
                push!(args, LogLevel(run[name_str]))
            elseif name_str == "compression_alg"
                push!(args, compression_alg_dict[run[name_str]])
            else
                push!(args, run[name_str])
            end
        end
        
        new(args...)
    end
end

#returns a dictionary that maps a string version of each enum option to its value
function get_enum_dict(enum_type::DataType)
    options = instances(enum_type)
    names = map(string, options)
    key_val_pairs = zip(names, options)

    Dict{String, enum_type}(key_val_pairs)
end

mutable struct Config
    run::Run
    rng::Random.MersenneTwister

    function Config(run::Run, seed_offset::UInt64)
        if run.fix_rng_seed
            seed = run.user_rng_seed + seed_offset
        else
            dev = Random.RandomDevice()
            seed = UInt64(Random.rand(dev) * 0xffffffffffffffff) + seed_offset
        end
            
        new(run, Random.MersenneTwister(seed))
    end
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
    
    Config(run, 0)
end

end
