module LcsMod

#use the Hunt-Szymanski algorithm
#https://en.wikipedia.org/wiki/Hunt%E2%80%93Szymanski_algorithm

struct Node
    i::Int64
    j::Int64
    next::Int64 #index of next node in link array
end

function lcs(a::String, b::String)
    lcs(Array{Char, 1}(a), Array{Char, 1}(b))
end

function diff(a::String, b::String)
    diff(Array{Char, 1}(a), Array{Char, 1}(b))
end

function diff(a::Array{T, 1}, b::Array{T, 1}) where T
    seq_info = lcs(a, b)

    #new
    added = Array{Int64, 1}()
    removed = Array{Int64, 1}()

    #v1
    # input_index = 1
    # @inbounds for (next_a, next_b) in seq_info
    #     next = min(next_a, next_b)
    #     #these items are absent from the lcs, but present in a, so they were removed
    #     @inbounds while input_index < length(a) || input_index < length(b)
    #         if input_index < next
    #             append!(removed, collect(input_index:next - 1))
    #             append!(added, collect(input_index:next - 1))
    #             input_index = next + 1 #skip over the item at index next
    #         else
    #             #note: nothing is appended if input_index >= length
    #             append!(removed, collect(input_index:length(a)))
    #             append!(added, collect(input_index:length(b)))
    #             input_index = length(a) #terminate the loop
    #         end
    #     end
    # end

    #v2
    # input_index = 1
    # @inbounds for (next_a, next_b) in seq_info
    #     next = min(next_a, next_b)
    #     indices = collect(input_index:next - 1)
    #     append!(removed, indices)
    #     append!(added, indices)
    #     input_index = next + 1 #skip over the item at index next
    # end

    # #note: nothing is appended if input_index >= length
    # if input_index < length(a) || input_index < length(b)
    #     append!(removed, collect(input_index:length(a)))
    #     append!(added, collect(input_index:length(b)))
    # end


    #old
    added = Array{Int64, 1}()
    removed = Array{Int64, 1}()

    a_index = 1
    b_index = 1
    @inbounds for (next_a, next_b) in seq_info
        #these items are absent from the lcs, but present in a, so they were removed
        @inbounds while a_index < min(next_a, next_b)
            push!(removed, a_index)
            a_index += 1
        end

        #these items are absent from the lcs, but present in b, so they were added
        @inbounds while b_index < min(next_a, next_b)
            push!(added, b_index)
            b_index += 1
        end

        a_index += 1
        b_index += 1
    end

    #any remaining elements of a are all absent from the lcs, but present in a, so they were removed
    @inbounds while a_index <= length(a)
        push!(removed, a_index)
        a_index += 1
    end

    #any remaining elements of b are all absent from the lcs, but present in b, so they were added
    @inbounds while b_index <= length(b)
        push!(added, b_index)
        b_index += 1
    end

    (added=added, removed=removed)
end

function lcs(a::Array{T, 1}, b::Array{T, 1}) where T
    #swap params so the shorter one (if any) is a
    swap = length(a) > length(b)
    if swap
        a, b = b, a
    end
    max_len = length(b)
    min_len = length(a)
    
    #note: thresh and link are 0-indexed in the paper, so we have to handle the 0-index edge cases
    thresh = fill(max_len + 1, min_len)
    link = Array{Union{Node, Nothing}, 1}(nothing, min_len)
    #link = Dict{Int64, Union{Node, Nothing}}()

    #build linked lists
    matchlist = init_matchlist(a, b)
    #println("matchlist: $(matchlist)")

    #init link
    #link[0] = nothing

    #compute successive thresh values
    @inbounds for i in 1:length(a)
        @inbounds for j in matchlist[i]
            #find k such that thresh[k-1] < j <= thresh[k] (use binary search)
            #note: after this, 1 <= k <= n
            # println("i: $i")
            # println("matchlines[i]: $(matchlist[i])")
            # println("j: $j")
            # println("thresh: $(thresh)")

            k = searchsortedfirst(thresh, j)
            #println("k: $(k)")
            
            #lower = k - 1 == 0 ? 0 : thresh[k - 1]
            #upper = thresh[k]
            #println("$(lower) < $(j) <= $(upper)")

            if j < thresh[k]
                thresh[k] = j
                link[k] = Node(i, j, k - 1)
            end

            #println("link: $(link)")
            #println()
        end
    end

    #recover longest common subsequence in reverse order
    #find largest k such that thresh[k] != n + 1
    k = length(thresh)
    @inbounds while k > 0 && thresh[k] == max_len + 1
        k -= 1
    end

    lcs_indices = Array{Tuple{Int64, Int64}, 1}()
    #ptr = link[k]
    ptr = (k == 0) ? nothing : link[k]
    @inbounds while ptr != nothing
        index = swap ? (ptr.j, ptr.i) : (ptr.i, ptr.j)
        push!(lcs_indices, index)
        #ptr = link[ptr.next]
        ptr = (ptr.next == 0) ? nothing : link[ptr.next]
    end

    #Iterators.Reverse(lcs_indices)
    reverse(lcs_indices)
end

function init_matchlist(a::Array{T, 1}, b::Array{T, 1}) where T
    #since there may be repeated elements, construct matchlist such that we can push either an array of indices, or a pointer (view / subarray) to an existing matchlist entry (list)
    # matchlist = Array{Union{Array{Int64, 1}, SubArray}, 1}()
    matchlist = Array{Array{Int64, 1}, 1}(undef, length(a))
    @inbounds for i in 1:length(a)
        # lookup = Dict{T, Int64}() #element of a -> existing index in matchlist
        # #find all indices (j's) in b where a[i] == b[j], in reverse order
        # if a[i] in keys(lookup)
        #     index = lookup[a[i]]
        #     push!(matchlist, view(matchlist, index)) #push a view (type is SubArray) to avoid the copy
        # else
        match_indices = similar(b, Int64)
        count = 1
        @inbounds for j in eachindex(b)
            match_indices[count] = j
            count += (a[i] == b[j])
        end

        #push!(matchlist, list)
        matchlist[i] = reverse(resize!(match_indices, count - 1))
            # lookup[a[i]] = length(matchlist)
        # end
    end

    matchlist
end

end
