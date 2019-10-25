module DataMod

using IndividualMod
using CellTreeMod
using RunMod
import Serialization
import CodecZlib

#default to reading the first run for now...
#later it might be useful to give the user a way to choose
function read_data()
    println("Reading data...")
    
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
    println("Done.")
    
    (run, ea_pops, reg_trees)
end

end
