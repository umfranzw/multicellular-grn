module CellTreeMod

using CellMod
using SymMod
using DataStructures: Queue, enqueue!, dequeue!, isempty

export traverse, find_empty, find

#deapth-first traversal
function traverse(f::Function, node::Cell)
    f(node)
    for child in node.children
        traverse(child, f)
    end
end

#breadth-first traversal
function traverse_bf(f::Function, node::Cell)
    q = Queue{Cell}()
    enqueue!(q, node)

    while !isempty(q)
        cell = dequeue!(q)
        f(cell)
        
        for child in cell.children
            enqueue!(q, child)
        end
    end
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

end
