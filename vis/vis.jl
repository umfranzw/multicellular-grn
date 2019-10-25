using Gtk
using Gadfly
using Printf
using RunMod

mutable struct ControlState
    indiv::Int64
    ea_step::Int64
    reg_step::Int64
end

function main()
    run, ea_pops, reg_trees = DataMod.read_data()
    
    win = GtkWindow("Vis", 400, 400)
    vbox = GtkBox(:v)
    push!(win, vbox)
    
    paned = GtkPaned(:h)
    push!(vbox, paned)

    push!(paned, build_regsim_pane())
    push!(paned, build_tree_pane())

    push!(vbox, build_control_area(run))

    condition = Condition()
    endit(w) = notify(condition)
    signal_connect(endit, win, :destroy)
    showall(win)
    wait(condition)
end

main()
