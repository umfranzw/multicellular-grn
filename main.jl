include("run_info.jl")

for run in RunInfo.get_config_channel()
    println(run)
end
