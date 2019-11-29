module DataMod

using IndividualMod
using CellTreeMod
using RunMod
import Serialization
import CodecZlib

#default to reading the first run for now...
#later it might be useful to give the user a way to choose
function read_data(filename::String)
    println("Reading data...")
    
    file_path = join((RunMod.DATA_PATH, filename), "/")

    #read in the data and deserialize it. This gives a 3-tuple of the form (run, ea_states, reg_states), where
    #to first item is a Run struct and the following two are dictionaries
    in_stream = open(file_path, "r")
    run, ea_states, reg_states = Serialization.deserialize(in_stream)
    close(in_stream)

    #decompress and deserialize the state dictionaries
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
