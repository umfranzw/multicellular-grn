using Gtk
using Gadfly
using Printf

struct ControlState
    indiv::Int64
    ea_step::Int64
    reg_step::Int64
end

function main()
    win = GtkWindow("Vis", 400, 400)
    vbox = GtkBox(:v)
    push!(win, vbox)
    
    paned = GtkPaned(:h)
    push!(vbox, paned)

    push!(paned, build_regsim_pane())
    push!(paned, build_tree_pane())

    push!(vbox, build_control_area())

    condition = Condition()
    endit(w) = notify(condition)
    signal_connect(endit, win, :destroy)
    showall(win)
    wait(condition)
end

function build_regsim_pane()
    pane = GtkBox(:v)
    push!(pane, GtkLabel("Reg Sim"))

    pane
end

function build_tree_pane()
    pane = GtkBox(:v)
    push!(pane, GtkLabel("Tree"))

    pane
end

function make_control(title::String, step_range::StepRange{Int64, Int64}, cs::ControlState, entry_sym::Symbol)
    grid = GtkGrid()
    
    label = GtkLabel(title)
    set_gtk_property!(label, :width_chars, 12)
    entry = GtkEntry()
    set_gtk_property!(entry, :text, string(step_range.start))
    set_gtk_property!(entry, :width_chars, 4)
    up_button = GtkButton("+")
    down_button = GtkButton("-")

    grid[1, 1] = label
    grid[2, 1] = entry
    grid[3, 1] = up_button
    grid[4, 1] = down_button

    signal_connect(button -> adjust_val(entry, step_range, 1, cs, entry_sym), up_button, :clicked)
    signal_connect(button -> adjust_val(entry, step_range, -1, cs, entry_sym), down_button, :clicked)

    grid
end

function adjust_val(entry::GtkEntry, step_range::StepRange{Int64, Int64}, dir::Int64, cs::ControlState, entry_sym::Symbol)
    cur = parse(Int64, get_gtk_property(entry, :text, String))
    cur = clamp(cur + dir * step_range.step, step_range.start, step_range.stop)
    set_gtk_property!(entry, :text, string(cur))

    
    
    update_graph()
end

function build_control_area()
    cs = ControlState(1, 1, 1)
    
    vbox = GtkBox(:v)
    push!(vbox, make_control("Individual: ", 1:2:10, cs))
    push!(vbox, make_control("EA Step: ", 1:2:10, cs))
    push!(vbox, make_control("Reg Step: ", 1:2:10, cs))

    vbox
end

function update_graph()
end

main()
