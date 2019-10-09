using RunMod
using EvAlgMod

function main()
    for run in RunInfo.get_config_channel()
        EvAlgMod.ev_alg(run)
    end
end

main()
