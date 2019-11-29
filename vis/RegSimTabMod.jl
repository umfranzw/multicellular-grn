module RegSimTabMod

using Gtk
using RunMod

mutable struct ControlState
    indiv::Int64
    ea_step::Int64
    reg_step::Int64
end

function build(run::Run)
    vbox = GtkBox(:v)
    
    paned = GtkPaned(:h)
    genome_pane = build_genome_pane()
    tree_pane = build_tree_pane()
    push!(paned, genome_pane)
    push!(paned, tree_pane)

    controls = build_control_area(run)

    push!(vbox, paned)
    push!(vbox, controls)

    vbox
end

function build_tree_pane()
    pane = GtkBox(:v)
    push!(pane, GtkLabel("Tree"))

    pane
end

function build_genome_pane()
    pane = GtkBox(:v)
    push!(pane, GtkLabel("Genome"))

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

    set_field!(cs, entry_sym, cur)
    
    update_graph(cs)
end

function build_control_area(run::Run)
    cs = ControlState(1, 1, 1)
    
    vbox = GtkBox(:v)
    indiv_range = 1:1:run.pop_size
    reg_range = 1:1:run.reg_steps
    push!(vbox, make_control("Individual: ", indiv_range, cs, :indiv))
    push!(vbox, make_control("EA Step: ", run.step_range, cs, :ea_step))
    push!(vbox, make_control("Reg Step: ", reg_range, cs, :reg_step))

    vbox
end

function update_graph(cs::ControlState)
end

end
