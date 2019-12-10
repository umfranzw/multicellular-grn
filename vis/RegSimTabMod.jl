module RegSimTabMod

using Gtk
using StatsPlots
using RunMod
using IndividualMod
using CellTreeMod
using ProteinStoreMod
using ProteinPropsMod
using TreeVisMod
import GeneMod
import MiscUtilsMod

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
    controls, controls_vbox, add_callbacks = build_control_area(run, reg_trees)

    #props area
    props_notebook = build_props_area(run, reg_trees, controls, add_callbacks)

    #tree graph
    tree_vbox = build_tree_graph_area(run, ea_pops, reg_trees, controls, add_callbacks)

    #genome_graph
    genome_vbox = build_genome_graph_area(run, ea_pops, reg_trees, controls, add_callbacks)

    #paned container for graphs
    paned = GtkPaned(:h)
    push!(paned, genome_vbox)
    push!(paned, tree_vbox)

    #hbox for control and props areas
    hbox = GtkBox(:h)
    push!(hbox, controls_vbox)
    push!(hbox, props_notebook)
    
    push!(outer_vbox, paned)
    push!(outer_vbox, hbox)

    outer_vbox
end

function build_control_area(
    run::Run,
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}}
)
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

    controls, controls_vbox, add_callbacks
end

function build_props_area(
    run::Run,
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    controls::Controls,
    add_callbacks::Dict{Symbol, Function}
)
    notebook = GtkNotebook()

    concs_tab = build_concs_tab(run, reg_trees, controls, add_callbacks)
    bindings_tab = build_bindings_tab(run, reg_trees, controls, add_callbacks)
    
    push!(notebook, concs_tab, "Protein Concs")
    push!(notebook, bindings_tab, "Gene Bindings")

    notebook
end

function build_concs_tab(
    run::Run,
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    controls::Controls,
    add_callbacks::Dict{Symbol, Function}
)
    conc_types = repeat([Float64], run.num_genes) #props, concs
    store = GtkListStore(String, conc_types...)

    populate_concs_store(store, reg_trees, controls)
    for (sym, connector) in add_callbacks
        connector(button -> populate_concs_store(store, reg_trees, controls))
    end

    props_view = GtkTreeView(GtkTreeModel(store))
    ren = GtkCellRendererText()
    props_col = GtkTreeViewColumn("Props", ren, Dict([("text", 0)]))
    push!(props_view, props_col)

    for i in 1:run.num_genes
        concs_col = GtkTreeViewColumn("Conc $i", ren, Dict([("text", i)]))
        push!(props_view, concs_col)
    end

    props_win = GtkScrolledWindow()
    set_gtk_property!(props_win, :hscrollbar_policy, Gtk.GConstants.GtkPolicyType.NEVER) #many Bothans died to bring us this information...
    push!(props_win, props_view)

    props_win
end

function build_bindings_tab(
    run::Run,
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    controls::Controls,
    add_callbacks::Dict{Symbol, Function}
)
    num_reg_cols = MiscUtilsMod.num_enum_vals(GeneMod.RegSites)
    num_prod_cols = MiscUtilsMod.num_enum_vals(GeneMod.ProdSites)
    store = GtkListStore(Int64, repeat([String], num_reg_cols + num_prod_cols)...) #gene index, reg site bindings, prod site bindings

    populate_bindings_store(store, reg_trees, controls)
    for (sym, connector) in add_callbacks
        connector(button -> populate_bindings_store(store, reg_trees, controls))
    end

    props_view = GtkTreeView(GtkTreeModel(store))
    ren = GtkCellRendererText()
    props_col = GtkTreeViewColumn("Gene", ren, Dict([("text", 0)]))
    push!(props_view, props_col)

    for i in 1:num_reg_cols
        enum_val = GeneMod.RegSites(i)
        label = split(MiscUtilsMod.enum_val_to_str(enum_val), ".")[2]
        concs_col = GtkTreeViewColumn(label, ren, Dict([("text", i)]))
        push!(props_view, concs_col)
    end

    for i in 1:num_prod_cols
        enum_val = GeneMod.ProdSites(i)
        label = split(MiscUtilsMod.enum_val_to_str(enum_val), ".")[2]
        concs_col = GtkTreeViewColumn(label, ren, Dict([("text", num_reg_cols + i)]))
        push!(props_view, concs_col)
    end

    props_win = GtkScrolledWindow()
    set_gtk_property!(props_win, :hscrollbar_policy, Gtk.GConstants.GtkPolicyType.NEVER)
    push!(props_win, props_view)

    props_win
end

function build_tree_graph_area(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    controls::Controls,
    add_callbacks::Dict{Symbol, Function}
)
    tree_vbox = GtkBox(:v)
    push!(tree_vbox, GtkLabel("Tree"))
    tree_graph = GtkImage()
    build_tree_plot(run, ea_pops, reg_trees, controls, tree_graph)
    push!(tree_vbox, tree_graph)
    show(tree_graph)

    for (sym, connector) in add_callbacks
        connector(button -> build_tree_plot(run, ea_pops, reg_trees, controls, tree_graph))
    end

    tree_vbox
end

function build_genome_graph_area(
    run::Run,
    ea_pops::Dict{String, Array{Array{Individual, 1}, 1}},
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    controls::Controls,
    add_callbacks::Dict{Symbol, Function}
)
    genome_vbox = GtkBox(:v)
    push!(genome_vbox, GtkLabel("Genome"))
    genome_graph = GtkImage()
    build_genome_plot(run, ea_pops, reg_trees, controls, genome_graph)
    push!(genome_vbox, genome_graph)
    show(genome_graph)

    for (sym, connector) in add_callbacks
        connector(button -> build_genome_plot(run, ea_pops, reg_trees, controls, genome_graph))
    end

    genome_vbox
end

function populate_bindings_store(
    store::GtkListStore,
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    controls::Controls
)
    while length(store) > 0
        pop!(store)
    end
    
    indiv_index, ea_step, reg_step, cell_index = map(sym -> get_control_val(getfield(controls, sym)), (:indiv, :ea_step, :reg_step, :cell))
    tree = reg_trees["after_reg_step"][ea_step][indiv_index][reg_step]
    cell = CellTreeMod.get_bf_node(tree, cell_index)

    for i in 1:length(cell.gene_states)
        row = Array{Any, 1}()
        push!(row, i)
        for site_type in instances(GeneMod.RegSites)
            site_index = Int64(site_type)
            bound_protein = cell.gene_states[i].reg_site_bindings[site_index]
            bind_str = bound_protein == nothing ? "-" : ProteinPropsMod.to_str(bound_protein.props)
            push!(row, bind_str)
        end

        for site_type in instances(GeneMod.ProdSites)
            site_index = Int64(site_type)
            prod_protein = cell.gene_states[i].prod_site_bindings[site_index]
            bind_str = prod_protein == nothing ? "-" : ProteinPropsMod.to_str(prod_protein.props)
            push!(row, bind_str)
        end

        push!(store, tuple(row...))
    end
end

function populate_concs_store(
    store::GtkListStore,
    reg_trees::Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}},
    controls::Controls
)
    while length(store) > 0
        pop!(store)
    end
    
    indiv_index, ea_step, reg_step, cell_index = map(sym -> get_control_val(getfield(controls, sym)), (:indiv, :ea_step, :reg_step, :cell))
    tree = reg_trees["after_reg_step"][ea_step][indiv_index][reg_step]
    cell = CellTreeMod.get_bf_node(tree, cell_index)
    proteins = ProteinStoreMod.get_all(cell.proteins)
    for protein in proteins
        props = ProteinPropsMod.to_str(protein.props)
        push!(store, (props, protein.concs...))
    end
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
        label = ProteinPropsMod.to_str(proteins[i].props)
        labels[:, i] = [label]
        concs[:, i] = proteins[i].concs
    end

    plot = groupedbar(concs, bar_position=:overlay, label=labels);
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
