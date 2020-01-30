module DataMod

using IndividualMod
using CellTreeMod
using RunMod
import TrackerMod
import Serialization
import CodecZlib

export Data

mutable struct Data
    #dictionary order is ea_step, index, reg_step
    trees::Dict{Tuple{Int64, Int64, Int64}, CellTree}
    #(ea_step, index, reg_step) => (position, size)
    trees_index::Dict{Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64}}
    #dictionary order is ea_step, index
    indivs::Dict{Tuple{Int64, Int64}, Individual}
    #(ea_step, index) => (position, size)
    indivs_index::Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64}}
    file_handle::IOStream
    run::Union{Run, Nothing}

    function Data(filename::String)
        file_path = join((RunMod.DATA_PATH, filename), "/")
        file_handle = open(file_path, "r")
        
        #ea_step, index, reg_step
        trees = Dict{Tuple{Int64, Int64, Int64}, CellTree}()
        trees_index = Dict{Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64}}()
        
        #ea_step, index
        indivs = Dict{Tuple{Int64, Int64}, Individual}()
        indivs_index = Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64}}()
        
        data = new(trees, trees_index, indivs, indivs_index, file_handle, nothing)
        data.run = get_run(data)
        create_index(data)

        data
    end
end

function close(data::Data)
    close(data.file_handle)
end

function get_run(data::Data)
    seek(data.file_handle, 0)
    tag_type = read(data.file_handle, TrackerMod.TagType)
    
    if tag_type != TrackerMod.RunState
        println("Error - the run is not the first thing in the data file.")
        close(data)
        exit(1)
    end
    
    size = read(data.file_handle, Int64)

    read_obj(data, position(data.file_handle), size)
end

function get_indiv(data::Data, ea_step::Int64, pop_index::Int64)
    key = (ea_step,  pop_index)
    get_obj(data, key, :indivs, :indivs_index)
end

function get_tree(data::Data, ea_step::Int64, pop_index::Int64, reg_step::Int64)
    key = (ea_step, pop_index, reg_step)
    get_obj(data, key, :trees, :trees_index)
end

function get_obj(data::Data, key::Any, store_field::Symbol, index_field::Symbol)
    store = getfield(data, store_field)
    index = getfield(data, index_field)
    
    if key ∉ keys(store)
        pos, size = index[key]
        store[key] = read_obj(data, pos, size)
    end    
    
    store[key]
end

function read_obj(data::Data, pos::Int64, size::Int64)
    seek(data.file_handle, pos)

    comp_obj = read(data.file_handle, size)
    ser_obj = CodecZlib.transcode(CodecZlib.GzipDecompressor, comp_obj)
    ser_buf = IOBuffer(ser_obj)

    Serialization.deserialize(ser_buf)
end

function create_index(data::Data)
    while !eof(data.file_handle)
        tag_type = read(data.file_handle, TrackerMod.TagType)
        if tag_type == TrackerMod.RunState
            #note: no tag here
            size = read(data.file_handle, Int64)
            #no need to add this to an index - it's always at the start of the data file
            
        elseif tag_type == TrackerMod.IndivState
            ea_step = read(data.file_handle, Int64)
            pop_index = read(data.file_handle, Int64)
            size = read(data.file_handle, Int64)
            
            key = (ea_step, pop_index)
            data.indivs_index[key] = (position(data.file_handle), size)

        elseif tag_type == TrackerMod.CellTreeState
            ea_step = read(data.file_handle, Int64)
            reg_step = read(data.file_handle, Int64)
            pop_index = read(data.file_handle, Int64)
            size = read(data.file_handle, Int64)

            key = (ea_step, pop_index, reg_step)
            data.trees_index[key] = (position(data.file_handle), size)
        end

        #go to next item
        seek(data.file_handle, position(data.file_handle) + size)
    end
end

end
