using Printf

function print_matrix(M)
    rows, cols = size(M)
    for i in 1:rows
        for j in 1:cols
            @printf("%5.2f", M[i, j])
        end
        @printf("\n")
    end
    @printf("\n")
end

function get(mat, row, col)
    height, width = size(mat)
    if row < 1 || col < 1 || row > height || col > width
        return 0.0
    else
        return mat[row, col]
    end
end

function test_2d()
    alpha = 0.1
    h = 0.75
    dt = 1
    rows = 9
    cols = 9
    timesteps = 5

    Tk = zeros((rows, cols))
    Tk[5, 5] = 1.0
    Tk1 = copy(Tk)

    print_matrix(Tk)

    for k in 1:timesteps
        for i in 1:rows
            for j in 1:cols
                Tk1[i, j] = (1 - 4 * dt * alpha / h^2) * Tk[i, j] + dt * alpha * ((get(Tk, i, j-1) + get(Tk, i-1, j) + get(Tk, i+1, j) + get(Tk, i, j+1)) / h^2)
            end
        end
        
        Tk = copy(Tk1)
        print_matrix(Tk)
    end
end

function test_1d()
    alpha = 0.1
    h = 0.75
    dt = 1
    cols = 9
    timesteps = 10

    Tk = zeros((1, cols))
    Tk[1, 5] = 1.0
    Tk1 = copy(Tk)

    print_matrix(Tk)

    for k in 1:timesteps
        for j in 1:cols
            Tk1[1, j] = (1 - 2 * dt * alpha / h^2) * Tk[1, j] + dt * alpha * ((get(Tk, 1, j-1) + get(Tk, 1, j + 1)) / h^2)
        end

        Tk = copy(Tk1)
        print_matrix(Tk)
    end
end

struct Node
    conc::Float64
    children::Array{Node, 1}
end

function build_tree()
    node = Node(0.5, [])

    nodes = map(copy, repeat(node, 7))

    root = nodes[1]
    l1l = nodes[2]
    l1r = nodes[3]
    push!(root.children, l1l, l1r)

    push!(l1l.children, nodes[4])
    push!(l1r.children, nodes[5:7])

    root
end

struct TreeInfo
    node_to_level::Dict{Node, Int64}
    level_to_node::Dict{Int64, Array{Node, 1}}

    function TreeInfo(root::Node)
        node_to_level = Dict{Node, Int64}()
        level_to_node = Dict{Int64, Array{Node, 1}}()
        if tree.root != nothing
            build_info(root, node_to_level, level_to_node, 0)
        end

        TreeInfo(node_to_level, level_to_node)
    end
end

function build_info(node::Node, node_to_level::Dict{Node, Int64}, level_to_node::Dict{Int64, Array{Node, 1}}, level::Int64)
    node_to_level[node] = level
    if level ∉ keys(level_to_node)
        level_to_node[level] = Array{Node}()
    end
    push!(level_to_node, node)

    for child in node.children
        build_info(child, node_to_level, level_to_node, level + 1)
    end
end

#note: this function will not handle diagonals (i.e. one of {level_offset, col_offset} must be set to 0)
function get_2D(node_level::Int64, node_col::Int64, level_offset::Int64, col_offset::Int64, info::TreeInfo)
    width = length(info.level_to_node[level])
    node = info.level_to_node[node_level][node_col]
    target_col = node_col + node_offset
    target_level = node_level + level_offset
    
    outside = target_col < 1 || target_col > width
    outside = outside || (level_offset == -1 && node.parent == nothing)
    outside = outside || (level_offset == 1 && length(node.children) == 0)

    if outside
        return 0.0
    else
        if row_offset == -1 && col_offset == 0
            num_siblings = length(node.parent.children)
            return node.parent.conc / num_siblings

        #note: there will be no edge connecting the nodes invoved here (they'll be siblings)
        elseif row_offset == 0 && col_offset != 0
            return info.level_to_node[target_level][target_col]
            
        elseif row_offset == 1 && col_offset == 0
            sum = 0.0
            for child in node.children
               sum += child.conc
            end
            return sum
        end
    end
end

#results dict is (level, col) -> new_conc
function diffuse_node(node::Node, info::TreeInfo, results::Dict{Tuple(Int64, Int64), Float64})
    diff_alpha = 0.1
    diff_h = 0.75
    diff_dt = 1.0
    
    level = info.node_to_level[node]
    col = indexof(node, info.level_to_node[level])[1]

    new_conc = (1 - 4 * diff_dt * diff_alpha / diff_h^2) *
        node.conc +
        diff_dt * diff_alpha * ((
            get_2D(level, col, 0, -1, info) +
            get_2D(level, col, -1, 0, info) +
            get_2D(level, col, 1, 0, info) +
            get_2D(level, col, 0, 1, info)
        ) / diff_h^2)

    results[(level, col)] = new_conc
end

#breadth-first traversal
function traverse_bf(f::Function, start_node::Node)
    q = Queue{Node}()
    enqueue!(q, start_node)

    while !isempty(q)
        node = dequeue!(q)
        f(cell)
        
        for child in node.children
            enqueue!(q, child)
        end
    end
end

function diffuse_tree(root::Node)
    info = TreeInfo(root)
    results = Dict{Tuple(Int64, Int64), Float64}()

    #apply the diffusion
    traverse_bf(node -> diffuse_node(node, info, results), root)

    #update the concs
    for ((level, col), new_conc) in results
        info.level_to_node[level][col].conc = new_conc
    end
end

function print_tree(root::Node)
    info = TreeInfo(root)
    for (level, nodes) in info.level_to_node
        for node in nodes
            @printf("%5.2f", node.conc)
        end
        print("\n")
    end
end

function test_tree()
    root = build_tree()
    print_tree(root)

    steps = 3
    
    for i in 1:steps
        diffuse_tree(root)
        print_tree(root)
    end
end

test_tree()
