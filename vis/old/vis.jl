using Gtk
using Printf
using RunMod
using DataMod
using RegSimTabMod
using ProteinAppActionsMod

function get_args()
    # if length(ARGS) != 1
    #     println("Usage: julia vis.jl <datafile>")
    #     exit(1)
    # end

    # ARGS[1]
    "data"
end

function main()
    ProteinAppActionsMod.init_app_actions()
    
    datafile = get_args()
    data = Data(datafile)
    
    win = GtkWindow("Vis", 1100, 600)
    vbox = GtkBox(:v)
    push!(win, vbox)
    
    reg_sim_tab = RegSimTabMod.build(data)
    push!(vbox, reg_sim_tab)

    end_condition = Condition()
    signal_connect(w -> DataMod.close(data), win, :destroy)
    signal_connect(w -> notify(end_condition), win, :destroy)
    showall(win)
    wait(end_condition)
end

main()
