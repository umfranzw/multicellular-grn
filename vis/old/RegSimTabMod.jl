module RegSimTabMod

using Gtk
using StatsPlots
using RunMod
using IndividualMod
using CellTreeMod
using ProteinStoreMod
using ProteinPropsMod
using TreeVisMod
using DataMod
using ChainMod
using ChainGraphMod
using CacheMod
import GeneMod
import MiscUtilsMod
import GtkUtilsMod

pixbuf_cache_size = 10

mutable struct Control
    title::String
    entry::GtkEntry
    range::StepRange{Int64, Int64}
end

mutable struct State
    indiv::Control
    ea_step::Control
    reg_step::Control
    cell::Control

    proteins_cache::Cache{Tuple{Int64, Int64, Int64, Int64}, GdkPixbuf}
    trees_cache::Cache{Tuple{Int64, Int64, Int64, Int64}, GdkPixbuf}
    chains_cache::Cache{Tuple{Int64, Int64, Int64, Int64}, GdkPixbuf}

    function State(indiv::Control, ea_step::Control, reg_step::Control, cell::Control)
        new(indiv, ea_step, reg_step, cell,
            Cache{Tuple{Int64, Int64, Int64, Int64}, GdkPixbuf}(RegSimTabMod.pixbuf_cache_size),
            Cache{Tuple{Int64, Int64, Int64, Int64}, GdkPixbuf}(RegSimTabMod.pixbuf_cache_size),
            Cache{Tuple{Int64, Int64, Int64, Int64}, GdkPixbuf}(RegSimTabMod.pixbuf_cache_size))
    end
end

function build(
    data::Data
)
    outer_vbox = GtkBox(:v)

    #control area
    state, state_vbox, add_callbacks = build_control_area(data)

    #props area
    props_notebook = build_props_area(data, state, add_callbacks)

    #graphs area
    graphs_notebook = build_graphs_area(data, state, add_callbacks)

    #hbox for control and props areas
    hbox = GtkBox(:h)
    push!(hbox, state_vbox)
    push!(hbox, props_notebook)
    
    push!(outer_vbox, graphs_notebook)
    push!(outer_vbox, hbox)

    outer_vbox
end

function build_control_area(
    data::Data
)
    state = State(
        Control("Individual: ", GtkEntry(), 1:1:data.run.pop_size), #indiv
        Control("EA Step: ", GtkEntry(), data.run.step_range), #ea
        Control("Reg Step: ", GtkEntry(), 1:1:data.run.reg_steps), #reg
        Control("Cell: ", GtkEntry(), 1:1:1) #cell
    )
    state_vbox = GtkBox(:v)
    
    control_syms = Array{Symbol, 1}([:indiv, :ea_step, :reg_step, :cell])
    add_callbacks = Dict{Symbol, Function}()
    for i in 1:length(control_syms)
        sym = control_syms[i]
        control = getfield(state, sym)
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

        up_callback = button -> adjust_control_val(control, 1, map(sym -> getfield(state, sym), control_syms[i + 1 : end]))
        down_callback = button -> adjust_control_val(control, -1, map(sym -> getfield(state, sym), control_syms[i + 1 : end]))
        update_cell_range_callback = button -> update_cell_range(state, data)
        
        signal_connect(up_callback, up_button, :clicked)
        signal_connect(update_cell_range_callback, up_button, :clicked)
        
        signal_connect(down_callback, down_button, :clicked)
        signal_connect(update_cell_range_callback, down_button, :clicked)

        add_callbacks[sym] = f -> map(button -> signal_connect(f, button, :clicked), (up_button, down_button))
        
        push!(state_vbox, grid)
    end
    update_cell_range(state, data)

    state, state_vbox, add_callbacks
end

function build_graphs_area(
    data::Data,
    state::State,
    add_callbacks::Dict{Symbol, Function}
)
    notebook = GtkNotebook()
    set_gtk_property!(notebook, :expand, true)

    protein_tab = build_protein_graph_area(data, state, add_callbacks)
    tree_tab = build_tree_graph_area(data, state, add_callbacks)
    chains_tab = build_chain_graph_area(data, state, add_callbacks)

    push!(notebook, protein_tab, "Protein Concentrations")
    push!(notebook, tree_tab, "Cell Tree")
    push!(notebook, chains_tab, "Interaction Chains")

    notebook
end

function build_props_area(
    data::Data,
    state::State,
    add_callbacks::Dict{Symbol, Function}
)
    notebook = GtkNotebook()
    set_gtk_property!(notebook, :expand, true)

    concs_tab = build_concs_tab(data, state, add_callbacks)
    initial_proteins_tab = build_initial_proteins_tab(data, state, add_callbacks)
    bindings_tab = build_bindings_tab(data, state, add_callbacks)
    
    push!(notebook, concs_tab, "Protein Concs")
    push!(notebook, initial_proteins_tab, "Initial Proteins")
    push!(notebook, bindings_tab, "Gene Bindings")

    notebook
end

function build_concs_tab(
    data::Data,
    state::State,
    add_callbacks::Dict{Symbol, Function}
)
    conc_types = repeat([Float64], data.run.num_genes) #props, concs
    store = GtkListStore(String, conc_types...)

    populate_concs_store(store, data, state)
    for (sym, connector) in add_callbacks
        connector(button -> populate_concs_store(store, data, state))
    end

    props_view = GtkTreeView(GtkTreeModel(store))
    ren = GtkCellRendererText()
    props_col = GtkTreeViewColumn("Props", ren, Dict([("text", 0)]))
    push!(props_view, props_col)

    for i in 1:data.run.num_genes
        concs_col = GtkTreeViewColumn("Conc $i", ren, Dict([("text", i)]))
        push!(props_view, concs_col)
    end

    props_win = GtkScrolledWindow()
    set_gtk_property!(props_win, :hscrollbar_policy, Gtk.GConstants.GtkPolicyType.NEVER) #many Bothans died to bring us this information...
    push!(props_win, props_view)

    props_win
end

function build_initial_proteins_tab(
    data::Data,
    state::State,
    add_callbacks::Dict{Symbol, Function}
)
    conc_types = repeat([Float64], data.run.num_genes) #props, concs
    store = GtkListStore(String, conc_types...)

    populate_initial_proteins_store(store, data, state)
    for sym in (:indiv, :ea_step)
        connector = add_callbacks[sym]
        connector(button -> populate_initial_proteins_store(store, data, state))
    end

    props_view = GtkTreeView(GtkTreeModel(store))
    ren = GtkCellRendererText()
    props_col = GtkTreeViewColumn("Props", ren, Dict([("text", 0)]))
    push!(props_view, props_col)

    for i in 1:data.run.num_genes
        concs_col = GtkTreeViewColumn("Conc $i", ren, Dict([("text", i)]))
        push!(props_view, concs_col)
    end

    props_win = GtkScrolledWindow()
    set_gtk_property!(props_win, :hscrollbar_policy, Gtk.GConstants.GtkPolicyType.NEVER) #many Bothans died to bring us this information...
    push!(props_win, props_view)

    props_win
end

function build_bindings_tab(
    data::Data,
    state::State,
    add_callbacks::Dict{Symbol, Function}
)
    num_reg_cols = length(instances(GeneMod.RegSite))
    num_prod_cols = length(instances(GeneMod.ProdSite))
    store = GtkListStore(Int64, repeat([String], num_reg_cols + num_prod_cols)...) #gene index, reg site bindings, prod site bindings

    populate_bindings_store(store, data, state)
    for (sym, connector) in add_callbacks
        connector(button -> populate_bindings_store(store, data, state))
    end

    props_view = GtkTreeView(GtkTreeModel(store))
    ren = GtkCellRendererText()
    props_col = GtkTreeViewColumn("Gene", ren, Dict([("text", 0)]))
    push!(props_view, props_col)

    for enum_val in instances(GeneMod.RegSite)
        int_val = Int64(enum_val)
        label = string(enum_val)
        concs_col = GtkTreeViewColumn(label, ren, Dict([("text", int_val)]))
        push!(props_view, concs_col)
    end

    for enum_val in instances(GeneMod.ProdSite)
        int_val = Int64(enum_val)
        label = string(enum_val)
        concs_col = GtkTreeViewColumn(label, ren, Dict([("text", num_reg_cols + int_val)]))
        push!(props_view, concs_col)
    end

    props_win = GtkScrolledWindow()
    set_gtk_property!(props_win, :hscrollbar_policy, Gtk.GConstants.GtkPolicyType.NEVER)
    push!(props_win, props_view)

    props_win
end

function build_tree_graph_area(
    data::Data,
    state::State,
    add_callbacks::Dict{Symbol, Function}
)
    tree_vbox = GtkBox(:v)
    push!(tree_vbox, GtkLabel("Tree"))
    tree_graph = GtkImage()
    build_tree_plot(data, state, tree_graph)
    push!(tree_vbox, tree_graph)

    for (sym, connector) in add_callbacks
        connector(button -> build_tree_plot(data, state, tree_graph))
    end

    win = GtkScrolledWindow()
    push!(win, tree_vbox)
    show(tree_graph)
    
    win
end

function build_protein_graph_area(
    data::Data,
    state::State,
    add_callbacks::Dict{Symbol, Function}
)
    protein_vbox = GtkBox(:v)
    push!(protein_vbox, GtkLabel("Protein Concentrations"))
    protein_graph = GtkImage()
    build_protein_plot(data, state, protein_graph)
    push!(protein_vbox, protein_graph)
    show(protein_graph)

    for (sym, connector) in add_callbacks
        connector(button -> build_protein_plot(data, state, protein_graph))
    end

    protein_vbox
end

function build_chain_graph_area(
    data::Data,
    state::State,
    add_callbacks::Dict{Symbol, Function}
)
    chain_vbox = GtkBox(:v)
    push!(chain_vbox, GtkLabel("Interaction Chains"))
    chain_graph = GtkImage()
    build_chain_plot(data, state, chain_graph)
    push!(chain_vbox, chain_graph)

    for (sym, connector) in add_callbacks
        connector(button -> build_chain_plot(data, state, chain_graph))
    end
    
    win = GtkScrolledWindow()
    push!(win, chain_vbox)
    show(chain_graph)

    win
end

function populate_bindings_store(
    store::GtkListStore,
    data::Data,
    state::State
)
    while length(store) > 0
        pop!(store)
    end
    
    indiv_index, ea_step, reg_step, cell_index = map(sym -> get_control_val(getfield(state, sym)), (:indiv, :ea_step, :reg_step, :cell))
    tree = DataMod.get_tree(data, ea_step, indiv_index, reg_step)
    cell = CellTreeMod.get_bf_node(tree, cell_index)

    for i in 1:length(cell.gene_states)
        row = Array{Any, 1}()
        push!(row, i)
        for site_type in instances(GeneMod.RegSite)
            site_index = Int64(site_type)
            bound_protein = cell.gene_states[i].reg_site_bindings[site_index]
            gene_seq = ProteinPropsMod.to_str(cell.gene_states[i].gene.reg_sites[site_index])
            protein_seq = bound_protein == nothing ? "-" : ProteinPropsMod.to_str(bound_protein.props)
            push!(row, "$(gene_seq) : $(protein_seq)")
        end

        for site_type in instances(GeneMod.ProdSite)
            site_index = Int64(site_type)
            bound_protein = cell.gene_states[i].prod_site_bindings[site_index]
            gene_seq = ProteinPropsMod.to_str(cell.gene_states[i].gene.prod_sites[site_index])
            protein_seq = bound_protein == nothing ? "-" : ProteinPropsMod.to_str(bound_protein.props)
            push!(row, "$(gene_seq) : $(protein_seq)")
        end

        push!(store, tuple(row...))
    end

    println("populate_bindings_store")
end

function populate_concs_store(
    store::GtkListStore,
    data::Data,
    state::State
)
    while length(store) > 0
        pop!(store)
    end
    
    indiv_index, ea_step, reg_step, cell_index = map(sym -> get_control_val(getfield(state, sym)), (:indiv, :ea_step, :reg_step, :cell))
    tree = DataMod.get_tree(data, ea_step, indiv_index, reg_step)
    cell = CellTreeMod.get_bf_node(tree, cell_index)
    proteins = ProteinStoreMod.get_all(cell.proteins)
    for protein in proteins
        props = ProteinPropsMod.to_str(protein.props)
        push!(store, (props, protein.concs...))
    end

    println("populate_concs_store")
end

function populate_initial_proteins_store(
    store::GtkListStore,
    data::Data,
    state::State
)
    while length(store) > 0
        pop!(store)
    end
    
    indiv_index, ea_step = map(sym -> get_control_val(getfield(state, sym)), (:indiv, :ea_step))
    indiv = DataMod.get_indiv(data, ea_step, indiv_index)
    for protein in indiv.initial_cell_proteins
        props = ProteinPropsMod.to_str(protein.props)
        push!(store, (props, protein.concs...))
    end

    println("populate_initial_proteins_store")
end

function build_tree_plot(
    data::Data,
    state::State,
    graph::GtkImage
)
    indiv_index, ea_step, reg_step, cell_index = map(sym -> get_control_val(getfield(state, sym)), (:indiv, :ea_step, :reg_step, :cell))

    key = (indiv_index, ea_step, reg_step, cell_index)
    if CacheMod.contains(state.trees_cache, key)
        pixbuf = state.trees_cache[key]
    else
        tree = DataMod.get_tree(data, ea_step, indiv_index, reg_step)
        pixbuf = TreeVisMod.plot(tree, cell_index)
        state.trees_cache[key] = pixbuf
    end
    set_gtk_property!(graph, :pixbuf, pixbuf)

    println("build_tree_plot")
end

function build_protein_plot(
    data::Data,
    state::State,
    graph::GtkImage
)
    #data: ea_step, indiv, reg_step
    indiv_index, ea_step, reg_step, cell_index = map(sym -> get_control_val(getfield(state, sym)), (:indiv, :ea_step, :reg_step, :cell))

    key = (indiv_index, ea_step, reg_step, cell_index)
    if CacheMod.contains(state.proteins_cache, key)
        pixbuf = state.proteins_cache[key]
    else
        tree = DataMod.get_tree(data, ea_step, indiv_index, reg_step)
        cell = CellTreeMod.get_bf_node(tree, cell_index)
        proteins = ProteinStoreMod.get_all(cell.proteins)

        concs = zeros(data.run.num_genes, length(proteins))
        labels = fill("", (1, length(proteins)))
        for i in 1:length(proteins)
            label = ProteinPropsMod.to_str(proteins[i].props)
            labels[:, i] = [label]
            concs[:, i] = proteins[i].concs
        end

        plot = groupedbar(concs, bar_position=:overlay, label=labels);
        pixbuf = GtkUtilsMod.plot_to_pixbuf(plot)
        state.proteins_cache[key] = pixbuf
    end
    set_gtk_property!(graph, :pixbuf, pixbuf)
    
    println("build_protein_plot")
end

function build_chain_plot(
    data::Data,
    state::State,
    graph::GtkImage
)
    indiv_index, ea_step, reg_step, cell_index = map(sym -> get_control_val(getfield(state, sym)), (:indiv, :ea_step, :reg_step, :cell))

    key = (indiv_index, ea_step, reg_step, cell_index)
    if CacheMod.contains(state.chains_cache, key)
        pixbuf = state.chains_cache[key]
    else
        chains = ChainMod.build_chain_graph(data, ea_step, indiv_index, cell_index)
        pixbuf = ChainGraphMod.plot(chains)
        state.chains_cache[key] = pixbuf
    end
    set_gtk_property!(graph, :pixbuf, pixbuf)

    println("build_chain_plot")
end

function update_cell_range(
    state::State,
    data::Data
)
    indiv_index, ea_step, reg_step, cell_index = map(sym -> get_control_val(getfield(state, sym)), (:indiv, :ea_step, :reg_step, :cell))
    tree = DataMod.get_tree(data, ea_step, indiv_index, reg_step)
    size = CellTreeMod.size(tree)
    state.cell.range = StepRange(state.cell.range.start, state.cell.range.step, size) # note: must construct a new range, since they're immutable
end

function adjust_control_val(
    control::Control,
    dir::Int64,
    reset::Array{Control, 1}
)
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
