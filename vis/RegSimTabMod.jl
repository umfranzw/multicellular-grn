module RegSimTabMod

using Gtk
using StatsPlots
using RunMod
using IndividualMod
using CellTreeMod
using ProteinStoreMod
using ProteinPropsMod
using TreeVisMod

genome_graph_path = "/tmp/genome_graph.png"
tree_graph_path = "/tmp/tree_graph.png"

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
end

mutable struct Graphs
    genome_graph::GtkImage
    tree_graph::GtkImage
end

function build(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}}
)
    vbox = GtkBox(:v)
    control_area, controls, graphs = build_control_area(run, ea_pops, reg_trees)
    props_area = build_props_area(run, ea_pops, reg_trees)
    
    paned = GtkPaned(:h)
    genome_pane = build_genome_pane(run, ea_pops, reg_trees, controls, graphs)
    tree_pane = build_tree_pane(run, ea_pops, reg_trees, controls, graphs)
    push!(paned, genome_pane)
    push!(paned, tree_pane)

    hbox = GtkBox(:h)
    push!(hbox, control_area)
    push!(hbox, props_area)
    
    push!(vbox, paned)
    push!(vbox, hbox)

    vbox
end

function build_tree_pane(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    controls::Controls,
    graphs::Graphs
)
    pane = GtkBox(:v)
    push!(pane, GtkLabel("Tree"))

    build_tree_plot(run, ea_pops, reg_trees, controls, graphs)
    push!(pane, graphs.tree_graph)
    show(graphs.tree_graph)

    pane
end

function build_genome_pane(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    controls::Controls,
    graphs::Graphs
)
    pane = GtkBox(:v)
    push!(pane, GtkLabel("Genome"))
    
    build_genome_plot(run, ea_pops, reg_trees, controls, graphs)
    push!(pane, graphs.genome_graph)
    show(graphs.genome_graph)
    
    pane
end

function build_tree_plot(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    controls::Controls,
    graphs::Graphs
)
    indiv_index, ea_step, reg_step, cell_index = map(sym -> get_control_val(getfield(controls, sym)), (:indiv, :ea_step, :reg_step, :cell))
    tree = reg_trees["after_reg_step"][ea_step][indiv_index][reg_step]

    TreeVisMod.gen_graph(tree, tree_graph_path, cell_index)
    set_gtk_property!(graphs.tree_graph, :file, tree_graph_path)
end

function build_genome_plot(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    controls::Controls,
    graphs::Graphs
)
    #reg_trees: ea_step, indiv, reg_step
    indiv_index, ea_step, reg_step, cell_index = map(sym -> get_control_val(getfield(controls, sym)), (:indiv, :ea_step, :reg_step, :cell))
    tree = reg_trees["after_reg_step"][ea_step][indiv_index][reg_step]
    cell = CellTreeMod.get_bf_node(tree, cell_index)
    proteins = ProteinStoreMod.get_all(cell.proteins)

    concs = zeros(run.num_genes, length(proteins))
    labels = fill("", (1, length(proteins)))
    for i in 1:length(proteins)
        buf = IOBuffer()
        ProteinPropsMod.show(buf, proteins[i].props)
        seek(buf, 0)
        label = chomp(read(buf, String)) #protein sequence string (remove the newline)

        labels[:, i] = [label]
        concs[:, i]= proteins[i].concs
    end

    plot = groupedbar(concs, bar_position=:stack, bar_width=0.5, label=labels);
    savefig(plot, genome_graph_path);
    set_gtk_property!(graphs.genome_graph, :file, genome_graph_path)
end

function setup_controls(
    run::Run,
    controls::Controls,
    graphs::Graphs,
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

        up_callback = button -> adjust_control_val(control, 1, map(sym -> getfield(controls, sym), control_syms[i + 1 : end]))
        down_callback = button -> adjust_control_val(control, -1, map(sym -> getfield(controls, sym), control_syms[i + 1 : end]))
        update_genome_graph_callback = button -> build_genome_plot(run, ea_pops, reg_trees, controls, graphs)
        update_tree_graph_callback = button -> build_tree_plot(run, ea_pops, reg_trees, controls, graphs)
        update_cell_range_callback = button -> update_cell_range(controls, reg_trees)

        signal_connect(up_callback, up_button, :clicked)
        signal_connect(update_genome_graph_callback, up_button, :clicked)
        signal_connect(update_tree_graph_callback, up_button, :clicked)
        signal_connect(update_cell_range_callback, up_button, :clicked)
        
        signal_connect(down_callback, down_button, :clicked)
        signal_connect(update_genome_graph_callback, down_button, :clicked)
        signal_connect(update_tree_graph_callback, down_button, :clicked)
        signal_connect(update_cell_range_callback, down_button, :clicked)

        push!(vbox, grid)
    end
    update_cell_range(controls, reg_trees)

    vbox
end

function update_cell_range(controls::Controls, reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}})
    indiv_index, ea_step, reg_step, cell_index = map(sym -> get_control_val(getfield(controls, sym)), (:indiv, :ea_step, :reg_step, :cell))
    tree = reg_trees["after_reg_step"][ea_step][indiv_index][reg_step]
    size = CellTreeMod.size(tree)
    controls.cell.range = StepRange(controls.cell.range.start, controls.cell.range.step, size) # note: must construct a new range, since they're immutable
end

function adjust_control_val(control::Control, dir::Int64, reset::Array{Control, 1})
    cur = parse(Int64, get_gtk_property(control.entry, :text, String))
    cur = clamp(cur + dir * control.range.step, control.range.start, control.range.stop)
    set_gtk_property!(control.entry, :text, string(cur))

    for control in reset
        reset_control_val(control)
    end
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
        Control("Cell: ", GtkEntry(), 1:1:1) #cell
    )
    graphs = Graphs(GtkImage(), GtkImage())
    vbox = setup_controls(run, controls, graphs, ea_pops, reg_trees)

    vbox, controls, graphs
end

function build_props_area(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}}
)
    vbox = GtkBox(:v)

    #row: props, concs
    conc_types = repeat([Float64], run.num_genes)
    store = GtkListStore(String, conc_types...)

    indiv_index, ea_step, reg_step, cell_index = map(sym -> get_control_val(getfield(controls, sym)), (:indiv, :ea_step, :reg_step, :cell))
    tree = reg_trees["after_reg_step"][ea_step][indiv_index][reg_step]
    cell = CellTreeMod.get_bf_node(tree, cell_index)
    proteins = ProteinStoreMod.get_all(cell.proteins)
    for protein in proteins
        buf = IOBuffer()
        show(buf, protein.props)
        seek(buf, 0)
        props = String(take!(buf))

        push!(store, (props, protein.concs...))
    end

    view = GtkTreeView(GtkTreeModel(store))
    
    ren = GtkCellRendererText()
    props_col = GtkTreeViewColumn("Props", ren, Dict([("text", 0)]))
    push!(view, props_col)

    for i in 1:run.num_genes
        concs_col = GtkTreeViewColumn("Conc $i", ren, Dict([("text", i)]))
        push!(view, concs_col)
    end

    push!(vbox, GtkLabel("Proteins"))
    push!(vbox, view)
    
    vbox
end

end
