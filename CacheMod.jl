module CacheMod

using DataStructures
import Base:getindex, setindex!, keys

export Cache

mutable struct Cache{K, V}
    size::Int64
    insert_order::Deque{K}
    data::Dict{K, V}

    function Cache{K, V}(size::Int64) where K where V
        new(size, Deque{K}(), Dict{K, V}())
    end
end

function getindex(cache::Cache{K, V}, key::K) where K where V
    getindex(cache.data, key)
end

function setindex!(cache::Cache{K, V}, val::V, key::K) where K where V
    #if we're not replacing an existing item (so the size will grow)
    if key ∉ keys(cache.data)
        if length(cache.data) + 1 > cache.size
            #remove the oldest item (key at back of insert_order)
            oldest_key = pop!(cache.insert_order)
            delete!(cache.data, oldest_key)
        end
        pushfirst!(cache.insert_order, key)
    end
    setindex!(cache.data, val, key)
end

keys(cache::Cache{K, V}) where K where V = keys(cache.data)

contains(cache::Cache{K, V}, key::K) where K where V = key in keys(cache.data)

end
