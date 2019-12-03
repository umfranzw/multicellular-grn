module RegSimTabMod

using Gtk
using RunMod
using DataFrames
import Gadfly
import Compose
using IndividualMod
using CellTreeMod
using ProteinStoreMod
using ProteinPropsMod

mutable struct ControlState
    indiv::Int64
    ea_step::Int64
    reg_step::Int64
    cell::Int64
end

function build(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}}
)
    vbox = GtkBox(:v)
    controls, control_state = build_control_area(run)
    
    paned = GtkPaned(:h)
    genome_pane = build_genome_pane(run, ea_pops, reg_trees, control_state)
    tree_pane = build_tree_pane(run, ea_pops, reg_trees, control_state)
    push!(paned, genome_pane)
    push!(paned, tree_pane)

    push!(vbox, paned)
    push!(vbox, controls)

    vbox
end

function build_tree_pane(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    cs::ControlState
)
    pane = GtkBox(:v)
    push!(pane, GtkLabel("Tree"))

    pane
end

function build_genome_pane(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    cs::ControlState
)
    pane = GtkBox(:v)
    push!(pane, GtkLabel("Genome"))

    plot = build_genome_plot(run, ea_pops, reg_trees, cs)
    graphic = Gadfly.render(plot)
    canvas = GtkCanvas(400, 300)
    push!(pane, canvas)
    show(canvas)

    Gtk.draw(canvas) do widget
        Compose.draw(Compose.CAIROSURFACE(canvas.back), graphic)
    end
    
    pane
end

function build_genome_plot(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    cs::ControlState
)
    Gadfly.set_default_plot_size(400, 300)

    #reg_trees: ea_step, indiv, reg_step
    tree = reg_trees["after_reg_step"][cs.ea_step][cs.indiv][cs.reg_step]
    cell = CellTreeMod.get_bf_node(tree, cs.cell)
    proteins = ProteinStoreMod.get_all(cell.proteins)

    labels = Array{String, 1}()
    concs = Array{Float64, 1}()
    for i in length(proteins)
        buf = IOBuffer()
        ProteinPropsMod.show(buf, proteins[i].props)
        seek(buf, 0)
        label = read(buf, String) #protein sequence string
        append!(labels, repeat([label], run.num_genes))
        append!(concs, proteins[i].concs)
    end
    cell_dframe = DataFrame(label=labels, concs=concs)

    plot = Gadfly.plot(
        Gadfly.Scale.y_continuous(minvalue=0, maxvalue=1),
        #proteins concs (lines)
        Gadfly.layer(
            cell_dframe,
            y=:concs,
            color=:label,
            Gadfly.Geom.bar(position=:stack)
        )
    )

    plot
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
    cs = ControlState(1, 1, 1, 1)
    
    vbox = GtkBox(:v)
    indiv_range = 1:1:run.pop_size
    reg_range = 1:1:run.reg_steps
    cell_range = 1:1:1
    push!(vbox, make_control("Individual: ", indiv_range, cs, :indiv))
    push!(vbox, make_control("EA Step: ", run.step_range, cs, :ea_step))
    push!(vbox, make_control("Reg Step: ", reg_range, cs, :reg_step))
    push!(vbox, make_control("Cell: ", cell_range, cs, :cell))

    vbox, cs
end

function update_graph(cs::ControlState)
end

end
