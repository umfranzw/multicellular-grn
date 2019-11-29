using Gtk
using Gadfly
using Printf
using RunMod
using DataMod
using RegSimTabMod

function main()
    run, ea_pops, reg_trees = DataMod.read_data(ARGS[1])
    
    win = GtkWindow("Vis", 400, 400)
    vbox = GtkBox(:v)
    push!(win, vbox)
    
    reg_sim_tab = RegSimTabMod.build(run)
    push!(vbox, reg_sim_tab)

    condition = Condition()
    endit(w) = notify(condition)
    signal_connect(endit, win, :destroy)
    showall(win)
    wait(condition)
end

main()
