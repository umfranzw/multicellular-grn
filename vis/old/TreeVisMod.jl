module TreeVisMod

using CellTreeMod
using CellMod
using DataStructures
using SymMod
using GraphVizMod

function gen_dot_code(tree::CellTree, cell_index::Int64)
    buf = IOBuffer()
    print(buf, "digraph G {\n")

    if tree.root != nothing
        #breadth-first traversal
        next_id = 1
        q = Queue{Tuple{Cell, Int64}}()
        enqueue!(q, (tree.root, next_id))
        next_id += 1

        while !isempty(q)
            cell, id = dequeue!(q)
            label = cell.sym == nothing ? "" : sym_str(cell.sym)
            print(buf, "$id [label=\"$label\"")
            if id == cell_index
                print(buf, "color=red, style=filled")
            end
            print(buf, "];\n")
            
            for child in cell.children
                enqueue!(q, (child, next_id))
                print(buf, "$id -> $next_id\n")
                next_id += 1
            end
        end
    end
    
    print(buf, "}")

    String(take!(buf))
end

function sym_str(sym::Sym)
    buf = IOBuffer()
    show(buf, sym)
    seek(buf, 0)
    
    String(take!(buf))
end

function plot(tree::CellTree, cell_index::Int64)
    dot_code = gen_dot_code(tree, cell_index)
    
    GraphVizMod.plot(dot_code)
end

end
