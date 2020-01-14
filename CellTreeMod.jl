module CellTreeMod

using CellMod
using SymMod
using DataStructures: Queue, enqueue!, dequeue!, isempty
using Printf
#import Base.deepcopy

export CellTree,
    traverse, find_empty, find, copy

mutable struct CellTree
    root::Union{Cell, Nothing}
end

# function deepcopy(tree::CellTree)
#     CellTree(copy(tree.root))
# end

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
    if tree.root != nothing
        count = size(tree.root, 0)
    end

    count
end

function size(node::Cell, count::Int64)
    count += 1

    for child in node.children
        count += size(child, count)
    end
    
    count
end

function to_expr_str(root::Cell)
    to_expr_str(root, "")
end

function to_expr_str(node::Cell, expr_str::String)
    sym = node.sym
    expr_str *= "$(sym.val)"
    if sym.type == SymMod.FcnCall
        expr_str *= "("
    end
    
    for i in 1:length(node.children)
        expr_str = to_expr_str(node.children[i], expr_str)
        if sym.type == SymMod.FcnCall && i < length(node.children)
            expr_str *= ", "
        else
            expr_str *= ")"
        end
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

function find(node::Cell; type::Union{SymType, Nothing}=nothing, val::Any=nothing)
    results = Array{Cell, 1}()

    found = node.sym != nothing
    if type != nothing
        found = found && node.sym.type == type
    end

    if val != nothing
        found = found && node.sym.val == val
    end
    
    if found
        push!(results, node)
    end

    for child in node.children
        push!(results, find(child, type=type, val=val)...)
    end

    results
end

function show(io::IO, tree::CellTree, ilevel::Int64=0)
    iprintln(io, "Cell Tree:", ilevel)
    
    num_cells = size(tree)
    iprint(io, @sprintf("size: %d", size), ilevel + 1)
end

end
