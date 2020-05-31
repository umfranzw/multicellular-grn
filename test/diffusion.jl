using Printf
using DataStructures

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

mutable struct Node
    concs::Array{Float64, 1}
    children::Array{Node, 1}
    parent::Union{Node, Nothing}
end

function build_tree(num_concs::Int64)
    n0 = Node(repeat([1.0], num_concs), [], nothing)
    n1 = Node(zeros(num_concs), [], nothing)
    n2 = Node(zeros(num_concs), [], nothing)
    n3 = Node(zeros(num_concs), [], nothing)
    n4 = Node(zeros(num_concs), [], nothing)
    n5 = Node(zeros(num_concs), [], nothing)
    n6 = Node(zeros(num_concs), [], nothing)

    push!(n0.children, n1, n2)
    n1.parent = n0
    n2.parent = n0

    push!(n1.children, n3)
    n3.parent = n1

    push!(n2.children, n4, n5, n6)
    n4.parent = n2
    n5.parent = n2
    n6.parent = n2

    n0
end

struct TreeInfo
    node_to_level::Dict{Node, Int64}
    level_to_node::Dict{Int64, Array{Node, 1}}

    function TreeInfo(root::Node)
        node_to_level = Dict{Node, Int64}()
        level_to_node = Dict{Int64, Array{Node, 1}}()
        if root != nothing
            build_info(root, node_to_level, level_to_node, 1)
        end

        new(node_to_level, level_to_node)
    end
end

function build_info(node::Node, node_to_level::Dict{Node, Int64}, level_to_node::Dict{Int64, Array{Node, 1}}, level::Int64)
    node_to_level[node] = level
    if level ∉ keys(level_to_node)
        level_to_node[level] = Array{Node, 1}()
    end
    push!(level_to_node[level], node)

    for child in node.children
        build_info(child, node_to_level, level_to_node, level + 1)
    end
end

#note: this function will not handle diagonals (i.e. one of {level_offset, col_offset} must be set to 0)
function get_2D(node_level::Int64, node_col::Int64, conc_col::Int64, row_offset::Int64, col_offset::Int64, info::TreeInfo, num_conc_cols::Int64)
    cur_node = info.level_to_node[node_level][node_col]
    result = 0.0

    #get from above
    if row_offset == -1 && col_offset == 0
        if cur_node.parent != nothing
            target_conc_col = conc_col + col_offset
            if target_conc_col != 0 && target_conc_col != num_conc_cols + 1 #if not falling off left or right edge
                num_siblings = length(cur_node.parent.children)
                result = cur_node.parent.concs[target_conc_col] / num_siblings
            end
        end

    #get from left
    elseif row_offset == 0 && col_offset == -1
        target_conc_col = conc_col + col_offset
        if target_conc_col == 0 #fall off the left edge
            #check to see if we have a left sibling we can grab from
            if node_col > 1 #have left sibling
                sibling = info.level_to_node[node_level][node_col - 1]
                result = sibling.concs[num_conc_cols] #return the rightmost conc
            end
        else #within bounds in current node
            result = cur_node.concs[target_conc_col]
        end

    #get from right
    elseif row_offset == 0 && col_offset == 1
        target_conc_col = conc_col + col_offset
        if target_conc_col == num_conc_cols + 1 #fall off the right edge
            #check to see if we have a right sibling we can grab from
            if node_col < length(info.level_to_node[node_level]) #have right sibling
                sibling = info.level_to_node[node_level][node_col + 1]
                result = sibling.concs[1] #return the leftmost conc
            end
        else #within bounds
            result = cur_node.concs[target_conc_col]
        end

    #get from below
    elseif row_offset == 1 && col_offset == 0
        if length(cur_node.children) > 0
            target_conc_col = conc_col + col_offset
            if target_conc_col != 0 && target_conc_col != num_conc_cols + 1 #if not falling off left or right edge
                sum = 0.0
                for child in cur_node.children
                    sum += child.concs[target_conc_col]
                end
                result = sum
            end
        end

    else
        error("Invalid offsets")
    end

    result
end

#results dict is (level, col) -> new_conc
function diffuse_node(node::Node, info::TreeInfo, results::Dict{Tuple{Int64, Int64}, Array{Float64, 1}}, num_concs::Int64)
    diff_alpha = 0.1
    diff_h = 0.75
    diff_dt = 1.0
    
    node_level = info.node_to_level[node]
    node_col = findall(n -> n == node, info.level_to_node[node_level])[1]

    new_concs = Array{Float64, 1}()
    for i in 1:num_concs
        new_conc = (1 - 4 * diff_dt * diff_alpha / diff_h^2) *
            node.concs[i] +
            diff_dt * diff_alpha * ((
                get_2D(node_level, node_col, i, 0, -1, info, num_concs) +
                get_2D(node_level, node_col, i, -1, 0, info, num_concs) +
                get_2D(node_level, node_col, i, 1, 0, info, num_concs) +
                get_2D(node_level, node_col, i, 0, 1, info, num_concs)
            ) / diff_h^2)
        push!(new_concs, new_conc)
    end

    results[(node_level, node_col)] = new_concs
end

#breadth-first traversal
function traverse_bf(f::Function, start_node::Node)
    q = Queue{Node}()
    enqueue!(q, start_node)

    while !isempty(q)
        node = dequeue!(q)
        f(node)
        
        for child in node.children
            enqueue!(q, child)
        end
    end
end

function diffuse_tree(root::Node, num_concs::Int64)
    info = TreeInfo(root)
    results = Dict{Tuple{Int64, Int64}, Array{Float64, 1}}()

    #apply the diffusion
    traverse_bf(node -> diffuse_node(node, info, results, num_concs), root)

    #update the concs
    for ((level, col), new_concs) in results
        #println(new_concs)
        info.level_to_node[level][col].concs = new_concs
    end
end

function print_tree(root::Node)
    info = TreeInfo(root)
    for level in 1:length(info.level_to_node)
        for node in info.level_to_node[level]
            @printf("%s, ", node.concs)
        end
        print("\n")
    end
    print("\n")
end

function test_1d()
    alpha = 0.25
    h = 1.0
    dt = 1
    cols = 9
    timesteps = 10

    println("alpha: $(alpha), h: $(h)")

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

function test_tree()
    steps = 3
    num_concs = 5
    
    root = build_tree(num_concs)
    print_tree(root)
    
    for i in 1:steps
        diffuse_tree(root, num_concs)
        print_tree(root)
    end
end

test_1d()
