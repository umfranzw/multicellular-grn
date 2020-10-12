module CheckpointMod

using DeepDiffs

export ChangeInfo

const chunk_size = 256 #bytes

mutable struct ChangeInfo
    removed::Array{UnitRange{Int64}, 1}
    added::Array{Tuple{Int64, Array{UInt8, 1}}, 1}

    function ChangeInfo()
        new(Array{UnitRange{Int64}, 1}(), Array{Tuple{Int64, Array{UInt8, 1}}, 1}())
    end
end

function compress(checkpoint::Array{UInt8, 1}, target::Array{UInt8, 1})
    first = 1
    last = chunk_size
    info = ChangeInfo()
    
    while last <= length(checkpoint) && last <= length(target)
        chunk1 = checkpoint[first:last]
        chunk2 = target[first:last]
        #chunk1 = view(checkpoint, first:last)
        #chunk2 = view(target, first:last)
        result = deepdiff(chunk1, chunk2)
        
        update_info(info, first, result, chunk2)

        first = last + 1
        last = first + chunk_size - 1
    end

    #target is smaller
    if length(target) < length(checkpoint)
        #finish last (non-full) chunk of target, if any
        last = length(target)
        if last >= first
            chunk1 = checkpoint[first:last]
            chunk2 = target[first:last]
            #chunk1 = view(checkpoint, first:last)
            #chunk2 = view(target, first:last)
            result = deepdiff(chunk1, chunk2)

            update_info(info, first, result, chunk2)
            first = last + 1
        end
        
        #finish remainder of checkpoint, if any
        if first <= length(checkpoint)
            push!(info.removed, first:length(checkpoint))
        end
        
    #checkpoint is smaller
    else
        #finish last (non-full) chunk of checkpoint, if any
        last = length(checkpoint)
        if last >= first
            chunk1 = checkpoint[first:last]
            chunk2 = target[first:last]
            #chunk1 = view(checkpoint, first:last)
            #chunk2 = view(target, first:last)
            result = deepdiff(chunk1, chunk2)

            update_info(info, first, result, chunk2)
            first = last + 1
        end

        #finish remainder of target, if any
        if first <= length(target)
            push!(info.added, (first, target[first:end]))
        end
    end

    info
end

function decompress(checkpoint::Array{UInt8, 1}, info::ChangeInfo)
    result = deepcopy(checkpoint)

    #do removals (proceed in reverse order, so the current iteration does not affect parts of the list the next one will operate on (since removed is in sorted order by start index)
    for range in Iterators.reverse(info.removed)
        deleteat!(result, range)
    end

    #do additions - no need to reverse here
    for (start, data) in info.added
        i = 0
        while i < length(data)
            insert!(result, start + i, data[i + 1])
            i += 1
        end
    end
    
    result
end

function update_info(info::ChangeInfo, first::Int64, result::VectorDiff{Array{UInt8,1}, Array{UInt8,1}}, chunk2::Array{UInt8, 1})
    offset = first - 1
    start = nothing
    stop = nothing
    for index in result.removed
        if stop == nothing
            start = index
            stop = index
        elseif index == stop + 1
            stop = index
        else #index > stop, start new range
            push!(info.removed, start + offset : stop + offset)
            start = index
            stop = index
        end
    end
    #push last one, if necessary
    if stop != nothing
        push!(info.removed, start + offset : stop + offset)
    end

    start = nothing
    stop = nothing
    for index in result.added
        if stop == nothing
            start = index
            stop = index
        elseif index == stop + 1
            stop = index
        else #index > stop, start new range
            push!(info.added, (start + offset, chunk2[start:stop]))
            start = index
            stop = index
        end
    end
    #push last one, if necessary
    if stop != nothing
        push!(info.added, (start + offset, chunk2[start:stop]))
    end
end

end
