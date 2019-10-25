using Gtk
using Gadfly
using Printf
using RunMod
using Serialization
using IndividualMod
using CellTreeMod
import CodecZlib

mutable struct ControlState
    indiv::Int64
    ea_step::Int64
    reg_step::Int64
end

function main()
    println("Reading data...")
    run, ea_pops, reg_trees = read_data()
    println("Done.")
    
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

#default to reading the first run for now...
#later it might be useful to give the user a way to choose
function read_data()
    #get the run so we can determine the data file path
    run = RunMod.get_first_run()
    file_path = join((RunMod.DATA_PATH, run.data_output_file), "/")

    #read in the data and deserialize it. This gives a 2-tuple of the form (ea_states, reg_states), where
    #each item is a dictionary
    in_stream = open(file_path, "r")
    ea_states, reg_states = Serialization.deserialize(in_stream)
    close(in_stream)

    #decompress everything
    #ea_states:
    ea_pops = Dict{String, Array{Array{Individual, 1}, 1}}()
    for (label, pops) in ea_states
        #label is a string, and pops is an array of arrays of UInt8 (bytes)
        for comp_pop in pops
            decomp_pop = CodecZlib.transcode(CodecZlib.GzipDecompressor, comp_pop)
            buf = IOBuffer(decomp_pop; read=true)
            deserial_pop = Serialization.deserialize(buf)
            
            if label ∉ keys(ea_pops)
                ea_pops[label] = Array{Array{Individual, 1}, 1}()
            end
            push!(ea_pops[label], deserial_pop)
        end
    end

    #reg_states:
    reg_trees = Dict{String, Array{Array{Array{CellTree, 1}, 1}, 1}}()
    for (label, trees) in reg_states
        for comp_trees in trees
            decomp = CodecZlib.transcode(CodecZlib.GzipDecompressor, comp_trees)
            buf = IOBuffer(decomp; read=true)
            deserial = Serialization.deserialize(buf)

            if label ∉ keys(reg_trees)
                reg_trees[label] = Array{Array{Array{CellTree, 1}, 1}, 1}()
            end
            push!(reg_trees[label], deserial)
        end
    end

    (run, ea_pops, reg_trees)
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

main()
