using Gtk
using Printf
using RunMod
using DataMod
using RegSimTabMod

function get_args()
    # if length(ARGS) != 1
    #     println("Usage: julia vis.jl <datafile>")
    #     exit(1)
    # end

    # ARGS[1]
    "data"
end

function main()
    datafile = get_args()
    run, ea_pops, reg_trees = DataMod.read_data(datafile)
    
    win = GtkWindow("Vis", 1100, 600)
    vbox = GtkBox(:v)
    push!(win, vbox)
    
    reg_sim_tab = RegSimTabMod.build(run, ea_pops, reg_trees)
    push!(vbox, reg_sim_tab)

    condition = Condition()
    endit(w) = notify(condition)
    signal_connect(endit, win, :destroy)
    showall(win)
    wait(condition)
end

main()
