module DataMod

using IndividualMod
using CellTreeMod
using RunMod
import TrackerMod
import Serialization
import CodecZlib

export Data

mutable struct IndexEntry
    positions::Array{Int64, 1}
    loaded::Bool

    function IndexEntry()
        new(Array{Int64, 1}(), false)
    end
end

mutable struct Data
    trees::Dict{Int64, Dict{Int64, Dict{Int64, CellTree}}}
    indivs::Dict{Int64, Dict{Int64, Individual}}
    index::Dict{Int64, IndexEntry}
    file_handle::IOStream
    run::Union{Run, Nothing}

    function Data(filename::String)
        file_path = join((RunMod.DATA_PATH, filename), "/")
        file_handle = open(file_path, "r")
        
        #ea_step, index, reg_step
        trees = Dict{Int64, Dict{Int64, Dict{Int64, CellTree}}}()
        
        #ea_step, index
        indivs = Dict{Int64, Dict{Int64, Individual}}()
        
        #ea_step => IndexEntry([chunk_start_positions], loaded?)
        index = Dict{Int64, IndexEntry}()

        data = new(trees, indivs, index, file_handle, nothing)
        create_index(data)

        data
    end
end

function close(data::Data)
    close(data.file_handle)
end

function create_index(data::Data)
    seek(data.file_handle, 0)
    while !eof(data.file_handle)
        chunk_start = position(data.file_handle)
        chunk_min_ea_step, chunk_max_ea_step, comp_chunk_size = read_next_header(data)
        for step in chunk_min_ea_step : chunk_max_ea_step
            if step != 0 #initial save is given ea_step of -1. The ea_steps start at 1. There is no iteration 0, but min_step can be set to -1.
                if step ∉ keys(data.index)
                    data.index[step] = IndexEntry()
                end
                push!(data.index[step].positions, chunk_start)
            end
        end

        #go to next chunk
        seek(data.file_handle, position(data.file_handle) + comp_chunk_size)
    end
end

function read_chunks_for_step(data::Data, ea_step::Int64)
    already_read = data.index[ea_step].loaded
    if !already_read
        positions = data.index[ea_step].positions
        for pos in positions
            seek(data.file_handle, pos)
            read_next_chunk(data)
        end
        #mark as loaded
        data.index[ea_step].loaded = true
    end
end

function read_all(data::Data)
    seek(data.file_handle, 0)
    while !eof(data.file_handle)
        read_next_chunk(data)
    end

    #mark all as loaded
    for ea_step in keys(data.index)
        data.index[ea_step].loaded = true
    end
end

function read_next_header(data::Data)
    chunk_min_ea_step = read(data.file_handle, Int64)
    chunk_max_ea_step = read(data.file_handle, Int64)
    comp_chunk_size = read(data.file_handle, Int64)

    return (chunk_min_ea_step, chunk_max_ea_step, comp_chunk_size)
end

function read_next_chunk(data::Data)
    chunk_min_ea_step, chunk_max_ea_step, comp_chunk_size = read_next_header(data)
    comp_chunk = read(data.file_handle, comp_chunk_size)
    decomp_chunk = CodecZlib.transcode(CodecZlib.GzipDecompressor, comp_chunk)

    chunk_buf = IOBuffer(decomp_chunk)
    while !eof(chunk_buf)
        tag_type = TrackerMod.TagType(read(chunk_buf, UInt8))
        serial_obj_size = read(chunk_buf, Int64)
        serial_obj_chunk = read(chunk_buf, serial_obj_size)

        obj_buf = IOBuffer(serial_obj_chunk)

        #this is always a tuple, the first element of which is a value from the TrackerMod.TagType enum
        obj = Serialization.deserialize(obj_buf)
        
        if tag_type == TrackerMod.RunState
            data.run = obj
        elseif tag_type == TrackerMod.RegState
            insert_tree(data, obj)
        elseif tag_type == TrackerMod.EAState
            insert_indiv(data, obj)
        end
    end
end

function insert_tree(data::Data, obj::Tuple{Int64, Int64, Int64, CellTree})
    ea_step, reg_step, index, tree = obj
    #dictionary order is ea_step, index, reg_step
    if ea_step ∉ keys(data.trees)
        data.trees[ea_step] = Dict{Int64, Dict{Int64, CellTree}}()
    end
    if index ∉ keys(data.trees[ea_step])
        data.trees[ea_step][index] = Dict{Int64, CellTree}()
    end
    data.trees[ea_step][index][reg_step] = tree
end

function insert_indiv(data::Data, obj::Tuple{Int64, Int64, Individual})
    ea_step, index, indiv = obj
    #dictionary order is ea_step, index
    if ea_step ∉ keys(data.indivs)
        data.indivs[ea_step] = Dict{Int64, Individual}()
    end
    data.indivs[ea_step][index] = indiv
end

end
