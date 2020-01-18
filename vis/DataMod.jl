module DataMod

using IndividualMod
using CellTreeMod
using RunMod
import TrackerMod
import Serialization
import CodecZlib

function read_data(filename::String)
    println("Reading data...")
    
    file_path = join((RunMod.DATA_PATH, filename), "/")

    #read in the data and deserialize it. This gives a 3-tuple of the form (run, ea_states, reg_states), where
    #to first item is a Run struct and the following two are dictionaries
    file_in = open(file_path, "r")

    run = nothing
    #ea_step, index, reg_step
    trees = Dict{Int64, Dict{Int64, Dict{Int64, CellTree}}}()
    #ea_step, index
    indivs = Dict{Int64, Dict{Int64, Individual}}()
    while !eof(file_in)
        comp_chunk_size = read(file_in, Int64)
        comp_chunk = read(file_in, comp_chunk_size)
        decomp_chunk = CodecZlib.transcode(CodecZlib.GzipDecompressor, comp_chunk)

        chunk_buf = IOBuffer(decomp_chunk)
        while !eof(chunk_buf)
            serial_obj_size = read(chunk_buf, Int64)
            serial_obj_chunk = read(chunk_buf, serial_obj_size)

            obj_buf = IOBuffer(serial_obj_chunk)

            #this is always a tuple, the first element of which is a value from the TrackerMod.TagType enum
            obj = Serialization.deserialize(obj_buf)
            tag_type = obj[1]
            #println(tag_type)
            
            if tag_type == TrackerMod.RunState
                run = obj[2]
            elseif tag_type == TrackerMod.RegState
                insert_tree(trees, obj)
            elseif tag_type == TrackerMod.EAState
                insert_indiv(indivs, obj)
            end
        end
    end
    
    close(file_in)

    run, indivs, trees
end

function insert_tree(trees::Dict{Int64, Dict{Int64, Dict{Int64, CellTree}}}, obj::Tuple{TrackerMod.TagType, Int64, Int64, Int64, CellTree})
    _, ea_step, reg_step, index, tree = obj
    #dictionary order is ea_step, index, reg_step
    if ea_step ∉ keys(trees)
        trees[ea_step] = Dict{Int64, Dict{Int64, CellTree}}()
    end
    if index ∉ keys(trees[ea_step])
        trees[ea_step][index] = Dict{Int64, CellTree}()
    end
    trees[ea_step][index][reg_step] = tree
end

function insert_indiv(indivs::Dict{Int64, Dict{Int64, Individual}}, obj::Tuple{TrackerMod.TagType, Int64, Int64, Individual})
    _, ea_step, index, indiv = obj
    #dictionary order is ea_step, index
    if ea_step ∉ keys(indivs)
        indivs[ea_step] = Dict{Int64, Individual}()
    end
    indivs[ea_step][index] = indiv
end

end
