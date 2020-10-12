using LcsMod
using DeepDiffs
import Blosc
import Random

function print_arrays(a::Array{UInt8, 1}, b::Array{UInt8, 1})
    i = 1
    while i < length(a) || i < length(b)
        if i < length(a) && i < length(b)
            println("$(i): $(a[i])\t$(b[i])")
        elseif i < length(a)
            println("$(i): $(a[i])")
        elseif i < length(b)
            println("$(i): \t\t$(b[i])")
        end
        i += 1
    end
end

function do_deepdiffs(a::Array{UInt8, 1}, b::Array{UInt8, 1}, chunk_size::Int64)
    added = Array{Int64, 1}()
    removed = Array{Int64, 1}()
    n = max(length(a), length(b))
    chunks = n ÷ chunk_size
    leftover = n % chunk_size
    for i in 1:chunks
        low = (i - 1) * chunk_size + 1
        high = i * chunk_size
        result = deepdiff(a[low:high], b[low:high])
        append!(removed, result.removed)
        append!(added, result.added)
    end
    if leftover > 0
        #note: in case where first element of range > second element, an empty array is returned - even if the range is outside the bounds of the array
        result = deepdiff(a[n - leftover + 1 : end], b[n - leftover + 1 : end])
        append!(removed, result.removed)
        append!(added, result.added)
    end

    added, removed

    #result = deepdiff(a, b)
    #result.added, result.removed
end

function do_lcs(a::Array{T, 1}, b::Array{T, 1}, chunk_size::Int64) where T
    added = Array{Int64, 1}()
    removed = Array{Int64, 1}()
    # local_a = reinterpret(UInt16, a)
    # local_b = reinterpret(UInt16, b)
    # chunk_size ÷= 2

    n = max(length(a), length(b))
    chunks = n ÷ chunk_size
    leftover = n % chunk_size
    for i in 1:chunks
        low = (i - 1) * chunk_size + 1
        high = i * chunk_size
        result = LcsMod.diff(a[low:high], b[low:high])
        append!(removed, result.removed)
        append!(added, result.added)
    end
    if leftover > 0
        #note: in case where first element of range > second element, an empty array is returned - even if the range is outside the bounds of the array
        result = LcsMod.diff(a[n - leftover + 1 : end], b[n - leftover + 1 : end])
        append!(removed, result.removed)
        append!(added, result.added)
    end

    added, removed
    
    #LcsMod.diff(a, b)
end

seed = 5
use_seed = true
dd_chunk_size = 256
dd_chunks = 44
lcs_chunk_size = 256
lcs_chunks = 44
if dd_chunks * dd_chunk_size != lcs_chunks * lcs_chunk_size
    throw(ErrorException("alg sizes not equal"))
end

n = lcs_chunks * lcs_chunk_size
num_changes = 20
println("n: $(n)")
println("num_changes: $(num_changes)")
println("lcs_chunk_size: $(lcs_chunk_size)")
println("lcs_chunks: $(lcs_chunks)")
println("dd_ chunk_size: $(dd_chunk_size)")
println("dd_chunks: $(dd_chunks)")
println()

if use_seed
    gen = Random.MersenneTwister(seed)
else
    gen = Random.MersenneTwister()
end

a = Random.rand(gen, UInt8, n)
b = deepcopy(a)
for i in 1:num_changes
    index = Random.rand(gen, 1:n)
    b[index] = Random.rand(gen, UInt8)
end
insert!(a, 10, 0x05)
insert!(a, 20, 0x05)

#a_big = collect(reinterpret(UInt16, a))
#b_big = collect(reinterpret(UInt16, b))

println("Timing LCS:")
@time lcs_added, lcs_removed = do_lcs(a, b, lcs_chunk_size)
println()

# comp_a = Blosc.compress(a)
# comp_b = Blosc.compress(b)
# big_comp_a = Blosc.compress(a_big)
# big_comp_b = Blosc.compress(b_big)
# println("comp_a: $(Base.summarysize(comp_a))")
# println("comp_b: $(Base.summarysize(comp_b))")
# println("big_comp_a: $(Base.summarysize(big_comp_a))")
# println("big_comp_b: $(Base.summarysize(big_comp_b))")

println("Timing DeepDiffs")
@time dd_added, dd_removed = do_deepdiffs(a, b, dd_chunk_size)
println()

println("lcs_removed:")
println(length(lcs_removed))
println("dd_removed:")
println(length(dd_removed))
println("lcs_added:")
println(length(lcs_added))
println("dd_added:")
println(length(dd_added))

println("lcs_removed:")
println(lcs_removed)
println("dd_removed:")
println(dd_removed)
println("lcs_added:")
println(lcs_added)
println("dd_added:")
println(dd_added)

