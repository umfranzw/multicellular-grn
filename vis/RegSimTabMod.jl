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

mutable struct Control
    title::String
    entry::GtkEntry
    range::StepRange{Int64, Int64}
end

mutable struct Controls
    indiv::Control
    ea_step::Control
    reg_step::Control
    cell::Control
    canvas::GtkCanvas
end

function build(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}}
)
    vbox = GtkBox(:v)
    control_area, controls = build_control_area(run, ea_pops, reg_trees)
    
    paned = GtkPaned(:h)
    genome_pane = build_genome_pane(run, ea_pops, reg_trees, controls)
    tree_pane = build_tree_pane(run, ea_pops, reg_trees, controls)
    push!(paned, genome_pane)
    push!(paned, tree_pane)

    push!(vbox, paned)
    push!(vbox, control_area)

    vbox
end

function build_tree_pane(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    controls::Controls
)
    pane = GtkBox(:v)
    push!(pane, GtkLabel("Tree"))

    pane
end

function build_genome_pane(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    controls::Controls
)
    pane = GtkBox(:v)
    push!(pane, GtkLabel("Genome"))

    plot = build_genome_plot(run, ea_pops, reg_trees, controls)
    graphic = Gadfly.render(plot)
    #canvas = GtkCanvas(400, 300)
    push!(pane, controls.canvas)
    show(controls.canvas)

    Gtk.draw(controls.canvas) do widget
        Compose.draw(Compose.CAIROSURFACE(controls.canvas.back), graphic)
    end
    
    pane
end

function build_genome_plot(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    controls::Controls
)
    Gadfly.set_default_plot_size(400, 300)

    #reg_trees: ea_step, indiv, reg_step
    indiv_index, ea_step, reg_step, cell_index = map(sym -> get_control_val(getfield(controls, sym)), (:indiv, :ea_step, :reg_step, :cell))
    tree = reg_trees["after_reg_step"][ea_step][indiv_index][reg_step]
    cell = CellTreeMod.get_bf_node(tree, cell_index)
    proteins = ProteinStoreMod.get_all(cell.proteins)
    println(proteins)

    labels = Array{String, 1}()
    concs = Array{Float64, 1}()
    genes = Array{Int64, 1}()
    for i in 1:length(proteins)
        buf = IOBuffer()
        ProteinPropsMod.show(buf, proteins[i].props)
        seek(buf, 0)
        label = read(buf, String) #protein sequence string
        append!(labels, repeat([label], run.num_genes))
        append!(concs, proteins[i].concs)
        append!(genes, 1:run.num_genes)
    end
    cell_dframe = DataFrame(label=labels, concs=concs, genes=genes)

    plot = Gadfly.plot(
        Gadfly.Scale.y_continuous(minvalue=0, maxvalue=1),
        #proteins concs (lines)
        Gadfly.layer(
            cell_dframe,
            x=:genes,
            y=:concs,
            color=:label,
            Gadfly.Geom.bar(position=:stack)
        )
    )

    plot
end

function setup_controls(
    controls::Controls,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}}
)
    vbox = GtkBox(:v)
    
    control_syms = Array{Symbol, 1}([:indiv, :ea_step, :reg_step, :cell])
    for i in 1:length(control_syms)
        control = getfield(controls, control_syms[i])
        grid = GtkGrid()
        
        label = GtkLabel(control.title)
        set_gtk_property!(label, :width_chars, 12)
        set_gtk_property!(control.entry, :text, string(control.range.start))
        set_gtk_property!(control.entry, :width_chars, 4)
        up_button = GtkButton("+")
        down_button = GtkButton("-")

        grid[1, 1] = label
        grid[2, 1] = control.entry
        grid[3, 1] = up_button
        grid[4, 1] = down_button

        up_callback = button -> adjust_control_val(control, 1, map(sym -> getfield(controls, sym), control_syms[i + 1 : end])) || update_genome_graph(controls, ea_pops, reg_trees)
        down_callback = button -> adjust_control_val(control, -1, map(sym -> getfield(controls, sym), control_syms[i + 1 : end])) || update_genome_graph(controls, ea_pops, reg_trees)
        
        signal_connect(up_callback, up_button, :clicked)
        signal_connect(down_callback, down_button, :clicked)

        push!(vbox, grid)
    end

    vbox
end

function adjust_control_val(control::Control, dir::Int64, reset::Array{Control, 1})
    cur = parse(Int64, get_gtk_property(control.entry, :text, String))
    cur = clamp(cur + dir * control.range.step, control.range.start, control.range.stop)
    set_gtk_property!(control.entry, :text, string(cur))

    for control in reset
        reset_control_val(control)
    end

    false
end

function reset_control_val(control::Control)
    set_gtk_property!(control.entry, :text, "1")
end

function get_control_val(control::Control)
    parse(Int64, get_gtk_property(control.entry, :text, String))
end

function build_control_area(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}}
)
    controls = Controls(
        Control("Individual: ", GtkEntry(), 1:1:run.pop_size), #indiv
        Control("EA Step: ", GtkEntry(), run.step_range), #ea
        Control("Reg Step: ", GtkEntry(), 1:1:run.reg_steps), #reg
        Control("Cell: ", GtkEntry(), 1:1:1), #cell
        GtkCanvas(400, 300)
    )
    vbox = setup_controls(controls, ea_pops, reg_trees)

    vbox, controls
end

function update_genome_graph(
    controls::Controls,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}}
)
    plot = build_genome_plot(run, ea_pops, reg_trees, controls)
    graphic = Gadfly.render(plot)
    show(controls.canvas)

    Gtk.draw(controls.canvas) do widget
       Compose.draw(Compose.CAIROSURFACE(controls.canvas.back), graphic)
    end

    # composition = Compose.compose(Compose.compose(Compose.context(), Compose.rectangle()), Compose.fill("tomato"))
    # Gtk.draw(controls.canvas) do widget
    #     Compose.draw(Compose.CAIROSURFACE(controls.canvas.back), composition)
    # end
end

end
