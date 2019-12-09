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

function build(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}}
)
    outer_vbox = GtkBox(:v)

    #control area
    controls = Controls(
        Control("Individual: ", GtkEntry(), 1:1:run.pop_size), #indiv
        Control("EA Step: ", GtkEntry(), run.step_range), #ea
        Control("Reg Step: ", GtkEntry(), 1:1:run.reg_steps), #reg
        Control("Cell: ", GtkEntry(), 1:1:1) #cell
    )
    controls_vbox = GtkBox(:v)
    
    control_syms = Array{Symbol, 1}([:indiv, :ea_step, :reg_step, :cell])
    add_callbacks = Dict{Symbol, Function}()
    for i in 1:length(control_syms)
        sym = control_syms[i]
        control = getfield(controls, sym)
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
        update_cell_range_callback = button -> update_cell_range(controls, reg_trees)
        
        signal_connect(up_callback, up_button, :clicked)
        signal_connect(update_cell_range_callback, up_button, :clicked)
        
        signal_connect(down_callback, down_button, :clicked)
        signal_connect(update_cell_range_callback, down_button, :clicked)

        add_callbacks[sym] = f -> map(button -> signal_connect(f, button, :clicked), (up_button, down_button))
        
        push!(controls_vbox, grid)
    end
    update_cell_range(controls, reg_trees)

    #props area
    props_vbox = GtkBox(:v)

    conc_types = repeat([Float64], run.num_genes) #props, concs
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

    push!(props_vbox, GtkLabel("Proteins"))
    push!(props_vbox, view)

    #tree graph
    tree_vbox = GtkBox(:v)
    push!(tree_vbox, GtkLabel("Tree"))
    tree_graph = GtkImage()
    build_tree_plot(run, ea_pops, reg_trees, controls, tree_graph)
    push!(tree_vbox, tree_graph)
    show(tree_graph)

    for (sym, connector) in add_callbacks
        connector(button -> build_tree_plot(run, ea_pops, reg_trees, controls, tree_graph))
    end

    #genome_graph
    genome_vbox = GtkBox(:v)
    push!(genome_vbox, GtkLabel("Genome"))
    genome_graph = GtkImage()
    build_genome_plot(run, ea_pops, reg_trees, controls, genome_graph)
    push!(genome_vbox, genome_graph)
    show(genome_graph)

    for (sym, connector) in add_callbacks
        connector(button -> build_genome_plot(run, ea_pops, reg_trees, controls, genome_graph))
    end

    #paned container for graphs
    paned = GtkPaned(:h)
    push!(paned, genome_vbox)
    push!(paned, tree_vbox)

    #hbox for control and props areas
    hbox = GtkBox(:h)
    push!(hbox, controls_vbox)
    push!(hbox, props_vbox)
    
    push!(outer_vbox, paned)
    push!(outer_vbox, hbox)

    outer_vbox
end

function build_tree_plot(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    controls::Controls,
    graph::GtkImage
)
    indiv_index, ea_step, reg_step, cell_index = map(sym -> get_control_val(getfield(controls, sym)), (:indiv, :ea_step, :reg_step, :cell))
    tree = reg_trees["after_reg_step"][ea_step][indiv_index][reg_step]

    TreeVisMod.gen_graph(tree, tree_graph_path, cell_index)
    set_gtk_property!(graph, :file, tree_graph_path)
end

function build_genome_plot(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    controls::Controls,
    graph::GtkImage
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
    set_gtk_property!(graph, :file, genome_graph_path)
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

end
