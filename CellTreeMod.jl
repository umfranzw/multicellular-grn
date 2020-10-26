module CellTreeMod

using CellMod
using SymMod
using DataStructures: Queue, enqueue!, dequeue!, isempty
using Printf
using MiscUtilsMod

import Base.show

export CellTree, TreeInfo,
    traverse, find_empty, find, copy

mutable struct CellTree
    root::Union{Cell, Nothing}

    function CellTree(root::Union{Cell, Nothing})
        new(root)
    end
end

struct TreeInfo
    cell_to_level::Dict{Cell, Int64}
    level_to_cell::Dict{Int64, Array{Cell, 1}}

    function TreeInfo(tree::CellTree)
        cell_to_level = Dict{Cell, Int64}()
        level_to_cell = Dict{Int64, Array{Cell, 1}}()
        if tree.root != nothing
            build_info(tree.root, cell_to_level, level_to_cell, 1)
        end

        new(cell_to_level, level_to_cell)
    end
end

function build_info(cell::Cell, cell_to_level::Dict{Cell, Int64}, level_to_cell::Dict{Int64, Array{Cell, 1}}, level::Int64)
    cell_to_level[cell] = level
    if level ∉ keys(level_to_cell)
        level_to_cell[level] = Array{Cell, 1}()
    end
    push!(level_to_cell[level], cell)

    for child in cell.children
        build_info(child, cell_to_level, level_to_cell, level + 1)
    end
end

function contains_syms(tree::CellTree, syms::Array{Sym, 1})
    result = false
    if tree.root != nothing
        result = contains_syms(tree.root, syms)
    end

    result
end

function contains_syms(cell::Cell, syms::Array{Sym, 1})
    result = cell.sym != nothing && cell.sym in syms
    i = 1
    while !result && i <= length(cell.children)
        result = contains_syms(cell.children[i], syms)
        i += 1
    end

    result
end

#depth-first traversal
function traverse(f::Function, start_node::Cell)
    f(start_node)
    for child in start_node.children
        traverse(f, child)
    end
end

function traverse(f::Function, tree::CellTree)
    if tree.root != nothing
        traverse(f, tree.root)
    end
end

#breadth-first traversal
function traverse_bf(f::Function, start_node::Cell)
    q = Queue{Cell}()
    enqueue!(q, start_node)

    while !isempty(q)
        cell = dequeue!(q)
        f(cell)
        
        for child in cell.children
            enqueue!(q, child)
        end
    end
end

function traverse_bf(f::Function, tree::CellTree)
    if tree.root != nothing
        traverse_bf(f, tree.root)
    end
end

function get_bf_node(tree::CellTree, index::Int64)
    if tree.root == nothing
        return nothing
    else
        return get_bf_node(tree.root, index)
    end
end

function get_bf_node(node::Cell, dist::Int64)
    if dist == 1
        return node
    else
        q = Queue{Cell}()
        for child in node.children
            enqueue!(q, child)
        end
        
        while !isempty(q) && dist > 1
            node = dequeue!(q)
            dist -= 1

            if dist > 1
                for child in node.children
                    enqueue!(q, child)
                end
            end
        end

        if dist == 1
            return node
        else
            return nothing
        end
    end
end

function size(tree::CellTree)
    count = 0
    CellTreeMod.traverse(cell -> count += 1, tree)

    count
end

function to_expr_str(root::Cell)
    to_expr_str(root, "")
end

function to_expr_str(node::Cell, expr_str::String)
    sym = node.sym
    if sym == nothing
        expr_str *= "_"
    else
        expr_str *= "$(sym.val)"
    end
    
    if sym != nothing && sym.type == SymMod.FcnCall
        expr_str *= "("
    end
    
    for i in 1:length(node.children)
        expr_str = to_expr_str(node.children[i], expr_str)
        if sym != nothing && (sym.type == SymMod.FcnCall && i < length(node.children))
            expr_str *= ", "
        end
    end
    
    if sym != nothing && sym.type == SymMod.FcnCall
        expr_str *= ")"
    end

    expr_str
end

function to_expr_str(tree::CellTree)
    if tree.root != nothing
        return to_expr_str(tree.root)
    end
end

function find_empty(node::Cell)
    results = Array{Cell, 1}()

    if node.sym == nothing
        push!(results, node)
    end

    for child in node.children
        push!(results, find_empty(child)...)
    end

    results
end

# function find(node::Cell; type::Union{SymType, Nothing}=nothing, val::Any=nothing)
#     results = Array{Cell, 1}()

#     found = node.sym != nothing
#     if type != nothing
#         found = found && node.sym.type == type
#     end

#     if val != nothing
#         found = found && node.sym.val == val
#     end

#     if found
#         push!(results, node)
#     end

#     for child in node.children
#         push!(results, find(child, type=type, val=val)...)
#     end

#     results
# end

function find_by_id(tree::CellTree, id::UInt64)
    result = nothing
    
    if tree.root != nothing
        result = find_by_id(tree.root, id)
    end

    return result
end

function find_by_id(cur::Cell, key::UInt64)
    result = nothing
    if cur.id == key
        result = cur
    else
        i = 1
        while result == nothing && i <= length(cur.children)
            result = find_by_id(cur.children[i], key)
            i += 1
        end
    end

    result
end

function show(io::IO, tree::CellTree, ilevel::Int64=0)
    iprintln(io, "CellTree", ilevel)
    if tree.root == nothing
        iprintln(io, "(nothing)", ilevel + 1)
    else
        show_tree(io, tree.root, ilevel + 1)
    end
end

function show_tree(io::IO, cell::Cell, ilevel::Int64)
    cell_sym_str = cell.sym == nothing ? "_" : string(cell.sym)[1:1] #note: the range keeps it a String (rather than a Char)
    cell_desc = "$(cell_sym_str) (age: $(cell.age))"
    iprint(io, cell_desc, ilevel)
    for child in cell.children
        show_tree(io, child, ilevel + 1)
    end
end

end
